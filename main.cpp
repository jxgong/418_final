// #include "common.h"
#include "physics.h"
#include "fileio.h"
#include <vector>
#include <math.h>
#include <float.h>
#include <iostream>

#include "CycleTimer.h"


// #define gravity 9.81f
// #define universal_gas_constant 287.05f
// #define boundary_idling_temp 375.f
// #define PHYSICS_H 1

#define first_deriv(varname,prev,curr,next,field,delta, default_value)\
                if (prev && next){ \
                    varname = ((next)->field - (prev)->field) / (2.f*(delta)); \
                } \
                else if (prev){ \
                    varname = (2.f*(default_value) - (curr).field - (prev)->field) / (2.f*(delta));\
                } \
                else if (next){ \
                    varname = ((next)->field + (curr).field - 2.f*(default_value)) / (2.f*(delta));\
                } \
                else{ \
                    varname = 0.f; \
                }

#define dim2Index(x,y,z,Lx,Ly) ((x) + ((y)*(Lx)) + ((z)*(Lx)*(Ly)))
#define nghbrsInd(x,y,z) ((x) + 3 * (y) + 9 * (z) + 13)

void simulateStepCuda(std::vector<Node>& new_nodes,
                      std::vector<Node>& nodes,
                      const stepParams params);
// StartupOptions parseOptions(int argc, char** argv);
// std::vector<Node> loadFromFile(std::string filename);
// bool saveToFile(std::vector<Node>, std::string filename);


float calculate_pressure(Node *node){
    float nV = 0.f;
    nV += node->rho_o2 / o2_molar_mass;
    nV += node->rho_n2 / n2_molar_mass;
    nV += node->rho_fuel / fuel_molar_mass;
    nV += node->rho_co2 / co2_molar_mass;
    nV += node->rho_nox / nox_molar_mass;
    nV += node->rho_h2o / h2o_molar_mass;
    return nV * 8.31 * node->temperature;
}

float temperature_to_viscosity(float temp){
    return 0.000001458f * sqrtf(temp*temp*temp) / (temp + 110.4f);
}

float internal_energy_to_temperature(float temperature){
    // TODO: enter empirical mapping of internal energy to temperature
    return boundary_idling_temp;
}

float o2_enthalpy(float temperature){
    float T = temperature;
    if (temperature < 700.f){
        return 1000.f*(31.332f*T-20.235f*T*T+57.866f*T*T*T-36.506f*T*T*T*T-0.007f/T-8.903f-0.f);
    }
    else if (temperature < 2000.f){
        return 1000.f*(30.032f*T+8.773f*T*T-3.998f*T*T*T+0.788f*T*T*T*T-0.742f/T-11.325f-0.f);
    }
    else{
        return 1000.f*(20.911f*T+10.721f*T*T-2.020f*T*T*T+0.146f*T*T*T*T+9.246f/T+5.338f-0.f);
    }
}

void simulateStep(std::vector<Node>& new_nodes,
                  std::vector<Node>& nodes,
                  const stepParams params){

    const int length = params.length;
    const int width = params.width;
    const int depth = params.depth;
    const float deltax = params.deltax;
    const float deltay = params.deltay;
    const float deltaz = params.deltaz;
    const float deltat = params.deltat;

    for (int k = 0; k < depth; k++){
        for (int j = 0; j < width; j++){
            for (int i = 0; i < length; i++){
                int index = i + length * (j + width * (k));
                if (!((nodes[index].temperature >= fuel_autoignition_point) || 
                      (nodes[index].temperature >= fuel_flash_point &&
                       params.sparks.begin() + index != params.sparks.end()))){ 
                           //TODO: check if ^this works
                    continue;
                }
                float rho_o2 = nodes[index].rho_o2;
                float nV_o2 = rho_o2 / o2_molar_mass;
                float rho_fuel = nodes[index].rho_fuel;
                float nV_fuel = rho_fuel / fuel_molar_mass;
                //TODO: Define nV_air

                float delta_o2, delta_n2, delta_fuel, delta_co2, delta_nox, delta_h2o;
                if (2.f*nV_o2 >= 25.f*nV_fuel){
                    // reaction is limited by fuel
                    delta_fuel = -reaction_rate_coefficient * nV_fuel * deltat;
                    delta_o2 = -reaction_rate_coefficient * nV_o2 * deltat;
                    delta_nox = -(2.f*delta_o2 - 25.f*delta_fuel);
                }
                else{
                    // reaction is limited by air
                    delta_o2 = -reaction_rate_coefficient * nV_o2 * deltat;
                    delta_fuel = delta_o2 / 12.5f;
                    delta_nox = 0.f;
                }
                delta_co2 = -delta_fuel * 8.f;
                delta_h2o = -delta_fuel * 9.f;
                delta_n2 = -delta_nox / 2.f;

                nodes[index].rho_o2 += delta_o2 * o2_molar_mass;
                nodes[index].rho_n2 += delta_n2 * n2_molar_mass;
                nodes[index].rho_fuel += delta_fuel * fuel_molar_mass;
                nodes[index].rho_co2 += delta_co2 * co2_molar_mass;
                nodes[index].rho_nox += delta_nox * nox_molar_mass;
                nodes[index].rho_h2o += delta_h2o * h2o_molar_mass;

                nodes[index].dQdt = -delta_o2 * o2_formation_enthalpy
                                    -delta_n2 * n2_formation_enthalpy
                                    -delta_fuel * fuel_formation_enthalpy
                                    -delta_co2 * co2_formation_enthalpy
                                    -delta_nox * nox_formation_enthalpy
                                    -delta_h2o * h2o_formation_enthalpy;

            }
        }
    }

    // case of second or later step
    for (int k = 0; k < depth; k++){
        for (int j = 0; j < width; j++){
            for (int i = 0; i < length; i++){
                int index = i + length * (j + width * (k));
                std::vector<Node*> nghbrs(27,NULL);
                const Node node = nodes[index];
                for (int n = 0; n < 27; n++){
                    int x,y,z;
                    x = (n % 3) - 1;
                    y = ((n/3) % 3) - 1;
                    z = (n/3)/3 - 1;
                    if (!((x<0 && i<=0)||(y<0 && j<=0)||(z < 0 && k<=0)||
                          (x>0 && i+1 >= length)||(y>0 && j+1 >= width)||(z>0 && k+1>=depth))){
                        nghbrs[n] = &nodes[index+x+y*length+z*length*width];
                    }
                }
                Node next_node;
                std::cout << node.rho_o2 << " " <<
                             node.rho_n2 << " " <<
                             node.rho_fuel << " " <<
                             node.rho_co2 << " " <<
                             node.rho_nox << " " << 
                             node.rho_h2o << " " <<
                             node.pressure << " " <<
                             node.temperature << " " <<
                             node.viscosity << " " <<
                             node.internal_energy << " " <<
                             node.u << " " <<
                             node.v << " " <<
                             node.w << std::endl;

                float rho = node.rho_o2+node.rho_n2+node.rho_fuel+node.rho_co2+node.rho_nox+node.rho_h2o;
                // std::cout << rho << " ";

                float dudx, dvdx, dwdx;
                float dudy, dvdy, dwdy;
                float dudz, dvdz, dwdz;
                first_deriv(dudx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],u,deltax,0.f);
                first_deriv(dvdx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],v,deltax,0.f);
                first_deriv(dwdx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],w,deltax,0.f);
                first_deriv(dudy,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],u,deltay,0.f);
                first_deriv(dvdy,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],v,deltay,0.f);
                first_deriv(dwdy,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],w,deltay,0.f);
                first_deriv(dudz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],u,deltaz,0.f);
                first_deriv(dvdz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],v,deltaz,0.f);
                first_deriv(dwdz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],w,deltaz,0.f);

                float drhodx, drhody, drhodz, drhodt;
                // update rho_o2
                first_deriv(drhodx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],rho_o2,deltax,0.f);
                first_deriv(drhody,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],rho_o2,deltay,0.f);
                first_deriv(drhodz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],rho_o2,deltaz,0.f);
                drhodt = -(rho*(dudx+dvdy+dwdz)+drhodx*node.u+drhody*node.v+drhodz*node.w);
                next_node.rho_o2 = node.rho_o2 + drhodt;
                // std::cout << drhodt << " ";
                // update rho_n2
                first_deriv(drhodx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],rho_n2,deltax,0.f);
                first_deriv(drhody,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],rho_n2,deltay,0.f);
                first_deriv(drhodz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],rho_n2,deltaz,0.f);
                drhodt = -(rho*(dudx+dvdy+dwdz)+drhodx*node.u+drhody*node.v+drhodz*node.w);
                next_node.rho_n2 = node.rho_n2 + drhodt;
                // std::cout << drhodt << " ";
                // update rho_fuel
                first_deriv(drhodx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],rho_fuel,deltax,0.f);
                first_deriv(drhody,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],rho_fuel,deltay,0.f);
                first_deriv(drhodz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],rho_fuel,deltaz,0.f);
                drhodt = -(rho*(dudx+dvdy+dwdz)+drhodx*node.u+drhody*node.v+drhodz*node.w);
                next_node.rho_fuel = node.rho_fuel + drhodt;
                // std::cout << drhodt << " ";
                // update rho_co2
                first_deriv(drhodx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],rho_co2,deltax,0.f);
                first_deriv(drhody,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],rho_co2,deltay,0.f);
                first_deriv(drhodz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],rho_co2,deltaz,0.f);
                drhodt = -(rho*(dudx+dvdy+dwdz)+drhodx*node.u+drhody*node.v+drhodz*node.w);
                next_node.rho_co2 = node.rho_co2 + drhodt;
                // std::cout << drhodt << " ";
                // update rho_nox
                first_deriv(drhodx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],rho_nox,deltax,0.f);
                first_deriv(drhody,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],rho_nox,deltay,0.f);
                first_deriv(drhodz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],rho_nox,deltaz,0.f);
                drhodt = -(rho*(dudx+dvdy+dwdz)+drhodx*node.u+drhody*node.v+drhodz*node.w);
                next_node.rho_nox = node.rho_nox + drhodt * deltat;
                // std::cout << drhodt << " ";
                // update rho_h2o
                first_deriv(drhodx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],rho_h2o,deltax,0.f);
                first_deriv(drhody,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],rho_h2o,deltay,0.f);
                first_deriv(drhodz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],rho_h2o,deltaz,0.f);
                drhodt = -(rho*(dudx+dvdy+dwdz)+drhodx*node.u+drhody*node.v+drhodz*node.w);
                next_node.rho_h2o = node.rho_h2o + drhodt * deltat;
                // std::cout << drhodt << std::endl;

                Node *r0c0 = nghbrs[nghbrsInd(-1,-1,0)];
                Node *r0c1 = nghbrs[nghbrsInd(-1,1,0)];
                Node *r1c0 = nghbrs[nghbrsInd(1,-1,0)];
                Node *r1c1 = nghbrs[nghbrsInd(1,1,0)];

                float d2udxdy = (r0c0 ? r0c0->u : 0.f)
                              - (r0c1 ? r0c1->u : 0.f)
                              - (r1c0 ? r1c0->u : 0.f)
                              + (r1c1 ? r1c1->u : 0.f);
                d2udxdy /= (4.f * deltax * deltay);
                float d2vdxdy = (r0c0 ? r0c0->v : 0.f)
                              - (r0c1 ? r0c1->v : 0.f)
                              - (r1c0 ? r1c0->v : 0.f)
                              + (r1c1 ? r1c1->v : 0.f);
                d2vdxdy /= (4.f * deltax * deltay);
                float d2wdxdy = (r0c0 ? r0c0->w : 0.f)
                              - (r0c1 ? r0c1->w : 0.f)
                              - (r1c0 ? r1c0->w : 0.f)
                              + (r1c1 ? r1c1->w : 0.f);
                d2wdxdy /= (4.f * deltax * deltay);

                r0c0 = nghbrs[nghbrsInd(-1,0,-1)];
                r0c1 = nghbrs[nghbrsInd(-1,0,1)];
                r1c0 = nghbrs[nghbrsInd(1,0,-1)];
                r1c1 = nghbrs[nghbrsInd(1,0,1)];
                float d2udxdz = (r0c0 ? r0c0->u : 0.f)
                              - (r0c1 ? r0c1->u : 0.f)
                              - (r1c0 ? r1c0->u : 0.f)
                              + (r1c1 ? r1c1->u : 0.f);
                d2udxdz /= (4.f * deltax * deltaz);
                float d2vdxdz = (r0c0 ? r0c0->v : 0.f)
                              - (r0c1 ? r0c1->v : 0.f)
                              - (r1c0 ? r1c0->v : 0.f)
                              + (r1c1 ? r1c1->v : 0.f);
                d2vdxdz /= (4.f * deltax * deltaz);
                float d2wdxdz = (r0c0 ? r0c0->w : 0.f)
                              - (r0c1 ? r0c1->w : 0.f)
                              - (r1c0 ? r1c0->w : 0.f)
                              + (r1c1 ? r1c1->w : 0.f);
                d2wdxdz /= (4.f * deltax * deltaz);

                r0c0 = nghbrs[nghbrsInd(0,-1,-1)];
                r0c1 = nghbrs[nghbrsInd(0,-1,1)];
                r1c0 = nghbrs[nghbrsInd(0,1,-1)];
                r1c1 = nghbrs[nghbrsInd(0,1,1)];
                float d2udydz = (r0c0 ? r0c0->u : 0.f)
                              - (r0c1 ? r0c1->u : 0.f)
                              - (r1c0 ? r1c0->u : 0.f)
                              + (r1c1 ? r1c1->u : 0.f);
                d2udydz /= (4.f * deltay * deltaz);
                float d2vdydz = (r0c0 ? r0c0->v : 0.f)
                              - (r0c1 ? r0c1->v : 0.f)
                              - (r1c0 ? r1c0->v : 0.f)
                              + (r1c1 ? r1c1->v : 0.f);
                d2vdydz /= (4.f * deltay * deltaz);
                float d2wdydz = (r0c0 ? r0c0->w : 0.f)
                              - (r0c1 ? r0c1->w : 0.f)
                              - (r1c0 ? r1c0->w : 0.f)
                              + (r1c1 ? r1c1->w : 0.f);
                d2wdydz /= (4.f * deltay * deltaz);

                Node *r0 = nghbrs[nghbrsInd(-1,0,0)];
                Node *r1 = nghbrs[nghbrsInd(0,0,0)];
                Node *r2 = nghbrs[nghbrsInd(1,0,0)];
                float d2udx2 = (r0 ? r0->u : 0.f) - 2.f * (r1 ? r1->u : 0.f) + (r2 ? r2->u : 0.f);
                d2udx2 /= (deltax * deltax);
                float d2vdx2 = (r0 ? r0->v : 0.f) - 2.f * (r1 ? r1->v : 0.f) + (r2 ? r2->v : 0.f);
                d2vdx2 /= (deltax * deltax);
                float d2wdx2 = (r0 ? r0->w : 0.f) - 2.f * (r1 ? r1->w : 0.f) + (r2 ? r2->w : 0.f);
                d2wdx2 /= (deltax * deltax);

                r0 = nghbrs[nghbrsInd(0,-1,0)];
                r1 = nghbrs[nghbrsInd(0,0,0)];
                r2 = nghbrs[nghbrsInd(0,1,0)];
                float d2udy2 = (r0 ? r0->u : 0.f) - 2.f * (r1 ? r1->u : 0.f) + (r2 ? r2->u : 0.f);
                d2udy2 /= (deltay * deltay);
                float d2vdy2 = (r0 ? r0->v : 0.f) - 2.f * (r1 ? r1->v : 0.f) + (r2 ? r2->v : 0.f);
                d2vdy2 /= (deltay * deltay);
                float d2wdy2 = (r0 ? r0->w : 0.f) - 2.f * (r1 ? r1->w : 0.f) + (r2 ? r2->w : 0.f);
                d2wdy2 /= (deltay * deltay);

                r0 = nghbrs[nghbrsInd(0,0,-1)];
                r1 = nghbrs[nghbrsInd(0,0,0)];
                r2 = nghbrs[nghbrsInd(0,0,1)];
                float d2udz2 = (r0 ? r0->u : 0.f) - 2.f * (r1 ? r1->u : 0.f) + (r2 ? r2->u : 0.f);
                d2udz2 /= (deltaz * deltaz);
                float d2vdz2 = (r0 ? r0->v : 0.f) - 2.f * (r1 ? r1->v : 0.f) + (r2 ? r2->v : 0.f);
                d2vdz2 /= (deltaz * deltaz);
                float d2wdz2 = (r0 ? r0->w : 0.f) - 2.f * (r1 ? r1->w : 0.f) + (r2 ? r2->w : 0.f);
                d2wdz2 /= (deltaz * deltaz);

                float dTxxdx = (2.f*d2udx2 - d2vdxdy - d2wdxdz) * 2.f / 3.f * node.viscosity;
                float dTyydy = (2.f*d2vdy2 - d2udxdy - d2wdydz) * 2.f / 3.f * node.viscosity;
                float dTzzdz = (2.f*d2wdz2 - d2udxdz - d2vdydz) * 2.f / 3.f * node.viscosity;
                float dTxydx = (d2vdx2 + d2udxdy) * node.viscosity;
                float dTxydy = (d2udy2 + d2vdxdy) * node.viscosity;
                float dTyzdy = (d2wdy2 + d2vdydz) * node.viscosity;
                float dTyzdz = (d2vdz2 + d2wdydz) * node.viscosity;
                float dTxzdx = (d2wdx2 + d2udxdz) * node.viscosity;
                float dTxzdz = (d2udz2 + d2wdxdz) * node.viscosity;

                float dpdx, dpdy, dpdz;
                first_deriv(dpdx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],pressure,deltax,node.pressure);
                first_deriv(dpdy,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],pressure,deltay,node.pressure);
                first_deriv(dpdz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],pressure,deltaz,node.pressure);

                float dudt = (-dpdx + dTxxdx + dTxydy + dTxzdz) / rho;
                float dvdt = (-dpdy + dTxydx + dTyydy + dTyzdz) / rho;
                float dwdt = (rho*gravity - dpdz + dTxzdx + dTyzdy + dTzzdz) / rho;

                // std::cout << "u " << node.u << " v " << node.v << " w " << node.w << " ";
                // std::cout << "dudt " << dudt << " dvdt " << dvdt << " dwdt " << dwdt << " ";
                next_node.u = node.u + dudt * deltat;
                next_node.v = node.v + dvdt * deltat;
                next_node.w = node.w + dwdt * deltat;
                // std::cout << "u " << next_node.u << " v " << next_node.v << " w " << next_node.w << std::endl;

                float dV = deltax*deltay*deltaz;
                float dEdx, dEdy, dEdz;
                first_deriv(dEdx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],internal_energy,deltax,boundary_specific_heat*boundary_idling_temp);
                first_deriv(dEdy,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],internal_energy,deltay,boundary_specific_heat*boundary_idling_temp);
                first_deriv(dEdz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],internal_energy,deltaz,boundary_specific_heat*boundary_idling_temp);
                dEdx *= dV;
                dEdy *= dV;
                dEdz *= dV;

                float Txx = (2.f*dudx - dvdy - dwdz) * 2.f/3.f *node.viscosity;
                float Tyy = (2.f*dvdy - dudx - dwdz) * 2.f/3.f *node.viscosity;
                float Tzz = (2.f*dwdz - dudx - dvdy) * 2.f/3.f *node.viscosity;
                float Txy = (dudy + dvdx) * node.viscosity;
                float Txz = (dudz + dwdx) * node.viscosity;
                float Tyz = (dvdz + dwdy) * node.viscosity;

                float dqxdx;
                float dkdx, dTdx, d2Tdx2;
                first_deriv(dkdx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],conductivity,deltax,boundary_conductivity);
                first_deriv(dTdx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],temperature,deltax,boundary_idling_temp);
                r0 = nghbrs[nghbrsInd(-1,0,0)];
                r1 = nghbrs[nghbrsInd(0,0,0)];
                r2 = nghbrs[nghbrsInd(1,0,0)];
                d2Tdx2 = (r0 ? r0->temperature : boundary_idling_temp) - 2.f * (r1 ? r1->u : boundary_idling_temp) + (r2 ? r2->u : boundary_idling_temp);
                d2Tdx2 /= (deltax * deltax);
                dqxdx = dkdx * dTdx + node.conductivity * d2Tdx2;

                float dqydy;
                float dkdy, dTdy, d2Tdy2;
                first_deriv(dkdy,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],conductivity,deltay,boundary_conductivity);
                first_deriv(dTdy,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],temperature,deltay,boundary_idling_temp);
                r0 = nghbrs[nghbrsInd(0,-1,0)];
                r1 = nghbrs[nghbrsInd(0,0,0)];
                r2 = nghbrs[nghbrsInd(0,1,0)];
                d2Tdy2 = (r0 ? r0->temperature : boundary_idling_temp) - 2.f * (r1 ? r1->u : boundary_idling_temp) + (r2 ? r2->u : boundary_idling_temp);
                d2Tdy2 /= (deltay * deltay);
                dqydy = dkdy * dTdy + node.conductivity * d2Tdy2;

                float dqzdz;
                float dkdz, dTdz, d2Tdz2;
                first_deriv(dkdz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],conductivity,deltax,boundary_conductivity);
                first_deriv(dTdz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],temperature,deltax,boundary_idling_temp);
                r0 = nghbrs[nghbrsInd(0,0,-1)];
                r1 = nghbrs[nghbrsInd(0,0,0)];
                r2 = nghbrs[nghbrsInd(0,0,1)];
                d2Tdz2 = (r0 ? r0->temperature : boundary_idling_temp) - 2.f * (r1 ? r1->u : boundary_idling_temp) + (r2 ? r2->u : boundary_idling_temp);
                d2udz2 /= (deltaz * deltaz);
                dqzdz = dkdz * dTdz + node.conductivity * d2Tdz2;

                float dEdt = node.dQdt;
                dEdt += rho*gravity*node.w;
                dEdt -= dEdx*node.u + dEdy*node.v + dEdz*node.w + 
                        node.internal_energy*dV*(dudx + dvdy + dwdz);
                dEdt -= dpdx*node.u + dpdy*node.v + dpdz*node.w +
                        node.pressure * (node.u+node.v+node.w);
                dEdt += dudx*Txx + dudy*Txy + dudz * Txz + 
                        node.u * (dTxxdx + dTxydy + dTxzdz);
                dEdt += dvdx*Txy + dvdy*Tyy + dvdz*Tyz +
                        node.v * (dTxydx + dTyydy + dTyzdz);
                dEdt += dwdx*Txz + dwdy*Tyz + dwdz*Tzz +
                        node.w * (dTxzdx + dTyzdy + dTzzdz);
                dEdt -= dqxdx + dqydy + dqzdz;

                next_node.internal_energy = node.internal_energy + dEdt / dV;
                next_node.temperature = internal_energy_to_temperature(next_node.internal_energy);
                next_node.viscosity = temperature_to_viscosity(next_node.temperature);
                next_node.conductivity = temperature_to_viscosity(next_node.temperature);
                next_node.pressure = calculate_pressure(&next_node);

                new_nodes[i] = next_node;
            }
            
        }
    }
    std::cout << std::endl;
    return;
}

int main(int argc, char *argv[]){
    StartupOptions options = parseOptions(argc, argv);
    stepParams params;
    int numIterations = 0;
    std::vector<Node> newNodes;
    std::vector<Node> nodes = loadFromFile(options.inputFile,&params,&numIterations);
    newNodes.resize(nodes.size());
    std::cout << params.width << " " << params.length << " " << params.depth << std::endl;
    for (int i = 0; i < numIterations; i++){
        double startTime = CycleTimer::currentSeconds();
        simulateStep(newNodes, nodes, params);
        double endTime = CycleTimer::currentSeconds();
        nodes.swap(newNodes);
    }
    saveToFile(newNodes, "output.txt");
    return 0;
}
