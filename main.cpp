#include "common.h"
#include "physics.h"
#include <vector>
#include <math.h>

#define first_deriv(varname,prev,curr,next,field,delta)\
                if (prev && next){ \
                    varname = (next->field - prev->field) / (2.f*delta); \
                } \
                else if (prev){ \
                    varname = (curr.field - prev->field) / (delta);\
                } \
                else if (next){ \
                    varname = (next->field - curr.field) / (delta);\
                } \
                else{ \
                    varneam = 0.f; \
                }

#define dim2Index(x,y,z,Lx,Ly) ((x) + ((y)*(Lx)) + ((z)*(Lx)*(Ly)))

float temperature(Node &node){
    return (node->rho_air+node->rho_fuel+node->rho_co2+node->rho_nox) 
            * 287.05f * node->pressure;
}

float viscocity(float temp){
    return 0.000001458f * sqrtf(temp*temp*temp) / (temp + 110.4f);
}

void simulateStep(std::vector<Node>& new_nodes,
                  const std::vector<Node>& nodes,
                  const stepParams){

    const int length = stepParams.length;
    const int width = stepParams.width;
    const int depth = stepParams.depth;
    const float deltax = stepParams.deltax;
    const float deltay = stepParams.deltay;
    const float deltaz = stepParams.deltaz;
    const float deltat = stepParams.deltat;

    // case of second or later step
    for (int k = 0; k < depth; k++){
        for (int j = 0; j < width; j++){
            for (int i = 0; i < length; i++){
                int index = i + length * (j + width * (k));
                const Node node = nodes[index];
                const Node *left = i > 0 ? &nodes[index-1] : NULL;
                const Node *right = i < length - 1 ? &nodes[index+1] : NULL;
                const Node *up = j > 0 ? &nodes[index-length] : NULL;
                const Node *down = j % width < width - 1 ? &nodes[index+length] : NULL;
                const Node *out = k > 0 ? &nodes[index-length*width] : NULL;
                const Node *in = k < depth -1 ? &nodes[index+length*width] : NULL;
                Node next_node;

                float dudx, dvdy,dwdz
                first_deriv(dudx,left,node,right,u,deltax);
                first_deriv(dvdy,up,node,down,v,deltay); = node.get_dudx(left,right,deltax);
                float dvdy = node.get_dudx(up,down,deltay);
                float dwdz = node.get_dwdz(in,out,deltaz);
                float drhodx, drhody, drhodz;
                // update rho_air
                first_deriv(drhodx,left,node,right,rho_air,deltax);
                first_deriv(drhody,up,node,down,rho_air,deltay);
                first_deriv(drhodz,out,node,in,rho_air,deltaz);
                drhodt = -(rho_curr*(dudx+dvdy+dwdz)+drhodx*(node.u+node.v+node.w));
                next_node.rho_air = node.rho_air + drhodt;
                // update rho_fuel
                first_deriv(drhodx,left,node,right,rho_fuel,deltax);
                first_deriv(drhody,up,node,down,rho_fuel,deltay);
                first_deriv(drhodz,out,node,in,rho_fuel,deltaz);
                drhodt = -(rho_curr*(dudx+dvdy+dwdz)+drhodx*(node.u+node.v+node.w));
                next_node.rho_fuel = node.rho_fuel + drhodt;
                // update rho_co2
                first_deriv(drhodx,left,node,right,rho_co2,deltax);
                first_deriv(drhody,up,node,down,rho_co2,deltay);
                first_deriv(drhodz,out,node,in,rho_co2,deltaz);
                drhodt = -(rho_curr*(dudx+dvdy+dwdz)+drhodx*(node.u+node.v+node.w));
                next_node.rho_co2 = node.rho_co2 + drhodt;
                // update rho_nox
                first_deriv(drhodx,left,node,right,rho_nox,deltax);
                first_deriv(drhody,up,node,down,rho_nox,deltay);
                first_deriv(drhodz,out,node,in,rho_nox,deltaz);
                drhodt = -(rho_curr*(dudx+dvdy+dwdz)+drhodx*(node.u+node.v+node.w));
                next_node.rho_nox = node.rho_nox + drhodt * deltat;

                int r0,r1,c0,c1;
                r0 = i > 0 ? i-1 : i;
                r1 = i < length-1 ? i+1 : i;
                c0 = j > 0 ? j-1 : j;
                c1 = j < width-1 ? j+1 : j;
                float d2udxdy = nodes[dim2Index(r0,c0,k,length,width)].u 
                              - nodes[dim2Index(r0,c1,k,length,width)].u
                              - nodes[dim2Index(r1,c0,k,length,width)].u
                              + nodes[dim2Index(r1,c1,k,length,width)].u;
                d2udxdy /= (((float) ((r1-r0)*(c1-c0))) * deltax * deltay);
                float d2vdxdy = nodes[dim2Index(r0,c0,k,length,width)].v 
                              - nodes[dim2Index(r0,c1,k,length,width)].v
                              - nodes[dim2Index(r1,c0,k,length,width)].v
                              + nodes[dim2Index(r1,c1,k,length,width)].v;
                d2vdxdy /= (((float) ((r1-r0)*(c1-c0))) * deltax * deltay);
                float d2wdxdy = nodes[dim2Index(r0,c0,k,length,width)].w 
                              - nodes[dim2Index(r0,c1,k,length,width)].w
                              - nodes[dim2Index(r1,c0,k,length,width)].w
                              + nodes[dim2Index(r1,c1,k,length,width)].w;
                d2wdxdy /= (((float) ((r1-r0)*(c1-c0))) * deltax * deltay);

                r0 = i > 0 ? i-1 : i;
                r1 = i < length-1 ? i+1 : i;
                c0 = k > 0 ? k-1 : k;
                c1 = k < depth-1 ? k+1 : k;
                float d2udxdz = nodes[dim2Index(r0,j,c0,length,width)].u 
                              - nodes[dim2Index(r0,j,c1,length,width)].u
                              - nodes[dim2Index(r1,j,c0,length,width)].u
                              + nodes[dim2Index(r1,j,c1,length,width)].u;
                d2udxdz /= (((float) ((r1-r0)*(c1-c0))) * deltax * deltaz);
                float d2vdxdz = nodes[dim2Index(r0,j,c0,length,width)].v 
                              - nodes[dim2Index(r0,j,c1,length,width)].v
                              - nodes[dim2Index(r1,j,c0,length,width)].v
                              + nodes[dim2Index(r1,j,c1,length,width)].v;
                d2vdxdz /= (((float) ((r1-r0)*(c1-c0))) * deltax * deltaz);
                float d2wdxdz = nodes[dim2Index(r0,j,c0,length,width)].w 
                              - nodes[dim2Index(r0,j,c1,length,width)].w
                              - nodes[dim2Index(r1,j,c0,length,width)].w
                              + nodes[dim2Index(r1,j,c1,length,width)].w;
                d2wdxdz /= (((float) ((r1-r0)*(c1-c0))) * deltax * deltaz);

                r0 = j > 0 ? j-1 : j;
                r1 = j < width-1 ? j+1 : j;
                c0 = k > 0 ? k-1 : k;
                c1 = k < depth-1 ? k+1 : k;
                float d2udydz = nodes[dim2Index(i,r0,c0,length,width)].u 
                              - nodes[dim2Index(i,r0,c1,length,width)].u
                              - nodes[dim2Index(i,r1,c0,length,width)].u
                              + nodes[dim2Index(i,r1,c1,length,width)].u;
                d2udxdz /= (((float) ((r1-r0)*(c1-c0))) * deltay * deltaz);
                float d2vdydz = nodes[dim2Index(i,r0,c0,length,width)].v 
                              - nodes[dim2Index(i,r0,c1,length,width)].v
                              - nodes[dim2Index(i,r1,c0,length,width)].v
                              + nodes[dim2Index(i,r1,c1,length,width)].v;
                d2vdxdz /= (((float) ((r1-r0)*(c1-c0))) * deltay * deltaz);
                float d2wdydz = nodes[dim2Index(i,r0,c0,length,width)].w 
                              - nodes[dim2Index(i,r0,c1,length,width)].w
                              - nodes[dim2Index(i,r1,c0,length,width)].w
                              + nodes[dim2Index(i,r1,c1,length,width)].w;
                d2wdydz /= (((float) ((r1-r0)*(c1-c0))) * deltay * deltaz);

                float d2udx2 = 0.f, d2vdx2 = 0.f , d2wdx2 = 0.f;
                if (i > 0 || i < length-1){
                    r0 = i == 0 ? i : (i == length-1 ? i-2 : i-1);
                    d2udx2 += nodes[dim2Index(r0,j,k,length,width)].u;
                    d2udx2 -= 2.f * nodes[dim2Index(r0+1,j,k,length,width)].u;
                    d2udx2 += nodes[dim2Index(r0+2,j,k,length,width)].u;
                    d2udx2 /= (deltax * deltax);
                    d2vdx2 += nodes[dim2Index(r0,j,k,length,width)].v;
                    d2vdx2 -= 2.f * nodes[dim2Index(r0+1,j,k,length,width)].v;
                    d2vdx2 += nodes[dim2Index(r0+2,j,k,length,width)].v;
                    d2vdx2 /= (deltax * deltax);
                    d2wdx2 += nodes[dim2Index(r0,j,k,length,width)].w;
                    d2wdx2 -= 2.f * nodes[dim2Index(r0+1,j,k,length,width)].w;
                    d2wdx2 += nodes[dim2Index(r0+2,j,k,length,width)].w;
                    d2wdx2 /= (deltax * deltax);
                }
                float d2udy2 = 0.f, d2vdy2 = 0.f, d2wdy2 = 0.f;
                if (j > 0 || j < width-1){
                    r0 = j == 0 ? j : (j == width-1 ? j-2 : j-1);
                    d2udy2 += nodes[dim2Index(i,r0,k,length,width)].u;
                    d2udy2 -= 2.f * nodes[dim2Index(i,r0+1,k,length,width)].u;
                    d2udy2 += nodes[dim2Index(i,r0+2,k,length,width)].u;
                    d2udy2 /= (deltay * deltay);
                    d2vdy2 += nodes[dim2Index(i,r0,k,length,width)].v;
                    d2vdy2 -= 2.f * nodes[dim2Index(i,r0+1,k,length,width)].v;
                    d2vdy2 += nodes[dim2Index(i,r0+2,k,length,width)].v;
                    d2vdy2 /= (deltay * deltay);
                    d2wdy2 += nodes[dim2Index(i,r0,k,length,width)].w;
                    d2wdy2 -= 2.f * nodes[dim2Index(i,r0+1,k,length,width)].w;
                    d2wdy2 += nodes[dim2Index(i,r0+2,k,length,width)].w;
                    d2wdy2 /= (deltay * deltay);
                }
                float d2udz2 = 0.f, d2vdz2 = 0.f, d2wdz2 = 0.f;
                if (k > 0 || k < depth-1){
                    r0 = k == 0 ? k : (k == depth-1 ? k-2 : k-1);
                    d2udz2 += nodes[dim2Index(i,j,r0,length,width)].u;
                    d2udz2 -= 2.f * nodes[dim2Index(i,j,r0+1,length,width)].u;
                    d2udz2 += nodes[dim2Index(i,j,r0+2,length,width)].u;
                    d2udz2 /= (deltaz * deltaz);
                    d2vdz2 += nodes[dim2Index(i,j,r0,length,width)].v;
                    d2vdz2 -= 2.f * nodes[dim2Index(i,j,r0+1,length,width)].v;
                    d2vdz2 += nodes[dim2Index(i,j,r0+2,length,width)].v;
                    d2vdz2 /= (deltaz * deltaz);
                    d2wdz2 += nodes[dim2Index(i,j,r0,length,width)].w;
                    d2wdz2 -= 2.f * nodes[dim2Index(i,j,r0+1,length,width)].w;
                    d2wdz2 += nodes[dim2Index(i,j,r0+2,length,width)].w;
                    d2wdz2 /= (deltaz * deltaz);
                }

                float dTxxdx = (2.f*d2udx2 - d2vdxdy - d2wdxdz) * 2.f / 3.f * node.viscocity;
                float dTyydy = (2.f*d2vdy2 - d2udxdy - d2wdydz) * 2.f / 3.f * node.viscocity;
                float dTzzdz = (2.f*d2wdz2 - d2udxdz - d2vdydz) * 2.f / 3.f * node.viscosity;
                float dTxydx = (d2vdx2 + d2udxdy) * node.viscosity;
                float dTxydy = (d2udy2 + d2vdxdy) * node.viscosity;
                float dTyzdy = (d2wdy2 + d2vdydz) * node.viscosity;
                float dTyzdz = (d2vdz2 + d2wdydz) * node.viscosity;
                float dTxzdx = (d2wdx2 + d2udxdz) * node.viscosity;
                float dTxzdz = (d2udz2 + d2wdxdz) * node.viscosity;

                float dpdx, dpdy, dpdz;
                first_deriv(dpdx,left,node,right,p,deltax);
                first_deriv(dpdy,up,node,down,p,deltay);
                first_deriv(dpdz,in,node,out,p,deltaz);

                float rho = node.rho_air+node.rho_fuel+node.rho_co2+node.rho_nox;

                float dudt = (-dpdx + dTxxdx + dTxydy + dTxzdz) / rho;
                float dvdt = (-dpdy + dTxydx + dTyydy + dTyzdz) / rho;
                float dwdt = (rho*g - dpdz + dTxzdx + dTyzdy + dTzzdz) / rho;

                next_node.u = node.u + dudt * deltat;
                next_node.v = node.v + dvdt * deltat;
                next_node.w = node.w + dwdt * deltat;

            }
        }
    }

    return;
}

int main(int argc, char *argv[]){
    StartupOptions options = parseOptions(argc, argv);
    return 0;
}