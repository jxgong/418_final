#include "common.h"
#include "physics.h"
#include <vector>
#include <math.h>

#define first_deriv(varname,prev,next,field,delta)\
                if (prev && next){ \
                    varname = ((next)->field - (prev)->field) / (2.f*(delta)); \
                } \
                else if (prev){ \
                    varname = -((prev)->field) / (2.f*(delta));\
                } \
                else if (next){ \
                    varname = ((next)->field) / (2.f*(delta));\
                } \
                else{ \
                    varname = 0.f; \
                }

#define dim2Index(x,y,z,Lx,Ly) ((x) + ((y)*(Lx)) + ((z)*(Lx)*(Ly)))
#define nghbrsInd(x,y,z) ((x) + 3 * (y) + 9 * (z))

float temperature(Node &node){
    return (node.rho_air+node.rho_fuel+node.rho_co2+node.rho_nox)
            * universal_gas_constant * node.pressure;
}

float viscosity(float temp){
    return 0.000001458f * sqrtf(temp*temp*temp) / (temp + 110.4f);
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

                float rho = node.rho_air+node.rho_fuel+node.rho_co2+node.rho_nox;
                float dQdt = 0.f;


                float dudx, dvdy, dwdz;
                first_deriv(dudx,nghbrs[nghbrsInd(-1,0,0)],nghbrs[nghbrsInd(1,0,0)],u,deltax);
                first_deriv(dvdy,nghbrs[nghbrsInd(0,-1,0)],nghbrs[nghbrsInd(0,1,0)],v,deltay);
                first_deriv(dwdz,nghbrs[nghbrsInd(0,0,-1)],nghbrs[nghbrsInd(0,0,1)],w,deltaz);

                float drhodx, drhody, drhodz, drhodt;
                // update rho_air
                first_deriv(drhodx,nghbrs[nghbrsInd(-1,0,0)],nghbrs[nghbrsInd(1,0,0)],rho_air,deltax);
                first_deriv(drhody,nghbrs[nghbrsInd(0,-1,0)],nghbrs[nghbrsInd(0,1,0)],rho_air,deltay);
                first_deriv(drhodz,nghbrs[nghbrsInd(0,0,-1)],nghbrs[nghbrsInd(0,0,1)],rho_air,deltaz);
                drhodt = -(rho*(dudx+dvdy+dwdz)+drhodx*node.u+drhody*node.v+drhodz*node.w);
                next_node.rho_air = node.rho_air + drhodt;
                // update rho_fuel
                first_deriv(drhodx,nghbrs[nghbrsInd(-1,0,0)],nghbrs[nghbrsInd(1,0,0)],rho_fuel,deltax);
                first_deriv(drhody,nghbrs[nghbrsInd(0,-1,0)],nghbrs[nghbrsInd(0,1,0)],rho_fuel,deltay);
                first_deriv(drhodz,nghbrs[nghbrsInd(0,0,-1)],nghbrs[nghbrsInd(0,0,1)],rho_fuel,deltaz);
                drhodt = -(rho*(dudx+dvdy+dwdz)+drhodx*node.u+drhody*node.v+drhodz*node.w);
                next_node.rho_fuel = node.rho_fuel + drhodt;
                // update rho_co2
                first_deriv(drhodx,nghbrs[nghbrsInd(-1,0,0)],nghbrs[nghbrsInd(1,0,0)],rho_co2,deltax);
                first_deriv(drhody,nghbrs[nghbrsInd(0,-1,0)],nghbrs[nghbrsInd(0,1,0)],rho_co2,deltay);
                first_deriv(drhodz,nghbrs[nghbrsInd(0,0,-1)],nghbrs[nghbrsInd(0,0,1)],rho_co2,deltaz);
                drhodt = -(rho*(dudx+dvdy+dwdz)+drhodx*node.u+drhody*node.v+drhodz*node.w);
                next_node.rho_co2 = node.rho_co2 + drhodt;
                // update rho_nox
                first_deriv(drhodx,nghbrs[nghbrsInd(-1,0,0)],nghbrs[nghbrsInd(1,0,0)],rho_nox,deltax);
                first_deriv(drhody,nghbrs[nghbrsInd(0,-1,0)],nghbrs[nghbrsInd(0,1,0)],rho_nox,deltay);
                first_deriv(drhodz,nghbrs[nghbrsInd(0,0,-1)],nghbrs[nghbrsInd(0,0,1)],rho_nox,deltaz);
                drhodt = -(rho*(dudx+dvdy+dwdz)+drhodx*node.u+drhody*node.v+drhodz*node.w);
                next_node.rho_nox = node.rho_nox + drhodt * deltat;

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
                first_deriv(dpdx,nghbrs[nghbrsInd(-1,0,0)],nghbrs[nghbrsInd(1,0,0)],pressure,deltax);
                first_deriv(dpdy,nghbrs[nghbrsInd(0,-1,0)],nghbrs[nghbrsInd(0,1,0)],pressure,deltay);
                first_deriv(dpdz,nghbrs[nghbrsInd(0,0,-1)],nghbrs[nghbrsInd(0,0,1)],pressure,deltaz);

                float dudt = (-dpdx + dTxxdx + dTxydy + dTxzdz) / rho;
                float dvdt = (-dpdy + dTxydx + dTyydy + dTyzdz) / rho;
                float dwdt = (rho*gravity - dpdz + dTxzdx + dTyzdy + dTzzdz) / rho;

                next_node.u = node.u + dudt * deltat;
                next_node.v = node.v + dvdt * deltat;
                next_node.w = node.w + dwdt * deltat;





                new_nodes[i] = next_node;
            }
        }
    }

    return;
}

int main(int argc, char *argv[]){
    StartupOptions options = parseOptions(argc, argv);
    return 0;
}
