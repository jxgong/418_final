#include "common.h"
#include "physics.h"
#include <vector>
#include "image.h"
#include "CycleTimer.h"
#include "cudaSimulator.h"
#include "ppm.h"
#include "fileio.h"

#include "math.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>


__constant__ int PRESSURE_LOWER_BOUND = 0.f;
__constant__ int PRESSURE_UPPER_BOUND = 10.f;
__constant__ int TEMP_UPPER_BOUND = 473.f;
__constant__ int TEMP_LOWER_BOUND = 273.f;
__constant__ int ENERGY_LOWER_BOUND = 0.f;
__constant__ int ENERGY_UPPER_BOUND = 1.f;
//TODO: check if these bounds are reasonable.

struct GlobalConstants{
    int imageWidth;
    int imageHeight;
    float* imageData;
    
    int nodeWidth, nodeLength, nodeDepth, sparkLen;
    float dx, dy, dz, dt;
    Node* nodes;
    Node* newNodes;
    Node* sparks;
};

__constant__ GlobalConstants cuImageData;

__device__ __inline__ int getNodeIdx(int x, int y, int z){
    return x + cuImageData.nodeLength * (y + cuImageData.nodeWidth * z);
}

__global__ void kernelSimStep(){
    int nodeX = blockIdx.x * blockDim.x + threadIdx.x;
    int nodeY = blockIdx.y * blockDim.y + threadIdx.y;
    int nodeZ = blockIdx.z * blockDim.z + threadIdx.z;
    int length = cuImageData.nodeLength;
    int width = cuImageData.nodeWidth;
    int depth = cuImageData.nodeDepth;
    if (nodeX >= length || nodeY >= width
                    || nodeZ >= depth){
        return;
    }
    float deltax = cuImageData.dx;
    float deltay = cuImageData.dy;
    float deltaz = cuImageData.dz;
    float deltat = cuImageData.dt;
    int index = getNodeIdx(nodeX, nodeY, nodeZ);
    const Node *node = &cuImageData.nodes[index];
    Node* nghbrs[27];
    for (int n = 0; n < 27; n++){
        int x,y,z;
        x = (n % 3) - 1;
        y = ((n/3) % 3) - 1;
        z = (n/9) - 1;
        if (!((x<0 && nodeX<=0)||(y<0 && nodeY<=0)||(z < 0 && nodeZ<=0)||
            (x>0 && nodeX+1 >= length)||(y>0 && nodeY+1 >= width)||
                (z>0 && nodeZ+1>=depth))){
            nghbrs[n] = &cuImageData.nodes[index+x+y*length+z*length*width];
        }
        else{
            nghbrs[n] = NULL;
        }
    }
    Node *next_node = &cuImageData.newNodes[index];
    float rho = node->rho_o2+node->rho_n2+node->rho_fuel+node->rho_co2+node->rho_nox+node->rho_h2o;
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

    float drhoudx, drhovdy, drhowdz, drhodt;
    // update rho_o2
    first_deriv_product(drhoudx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],rho_o2,u,deltax,0.f);
    first_deriv_product(drhovdy,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],rho_o2,v,deltay,0.f);
    first_deriv_product(drhowdz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],rho_o2,w,deltaz,0.f);
    drhodt = -drhoudx - drhovdy - drhowdz;//-(node->rho_o2*(dudx+dvdy+dwdz)+drhoudx*node->u+drhovdy*node->v+drhowdz*node->w);
    next_node->rho_o2 = node->rho_o2 + drhodt * deltat;
    if (next_node->rho_o2 < 0.f){
        next_node->rho_o2 = 0.f;
    }

    // update rho_n2
    first_deriv_product(drhoudx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],rho_n2,u,deltax,0.f);
    first_deriv_product(drhovdy,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],rho_n2,v,deltay,0.f);
    first_deriv_product(drhowdz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],rho_n2,w,deltaz,0.f);
    drhodt = -drhoudx - drhovdy - drhowdz;//drhodt = -(node->rho_n2*(dudx+dvdy+dwdz)+drhoudx*node->u+drhovdy*node->v+drhowdz*node->w);
    next_node->rho_n2 = node->rho_n2 + drhodt * deltat;
    if (next_node->rho_n2 < 0.f){
        next_node->rho_n2 = 0.f;
    }

    // update rho_fuel
    first_deriv_product(drhoudx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],rho_fuel,u,deltax,0.f);
    first_deriv_product(drhovdy,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],rho_fuel,v,deltay,0.f);
    first_deriv_product(drhowdz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],rho_fuel,w,deltaz,0.f);
    drhodt = -drhoudx - drhovdy - drhowdz;//drhodt = -(node->rho_fuel*(dudx+dvdy+dwdz)+drhoudx*node->u+drhovdy*node->v+drhowdz*node->w);
    next_node->rho_fuel = node->rho_fuel + drhodt * deltat;
    if (next_node->rho_fuel < 0.f){
        next_node->rho_fuel = 0.f;
    }

    // update rho_co2
    first_deriv_product(drhoudx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],rho_co2,u,deltax,0.f);
    first_deriv_product(drhovdy,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],rho_co2,v,deltay,0.f);
    first_deriv_product(drhowdz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],rho_co2,w,deltaz,0.f);
    drhodt = -drhoudx - drhovdy - drhowdz;//drhodt = -(node->rho_co2*(dudx+dvdy+dwdz)+drhoudx*node->u+drhovdy*node->v+drhowdz*node->w);
    next_node->rho_co2 = node->rho_co2 + drhodt * deltat;
    if (next_node->rho_co2 < 0.f){
        next_node->rho_co2 = 0.f;
    }

    // update rho_nox
    first_deriv_product(drhoudx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],rho_nox,u,deltax,0.f);
    first_deriv_product(drhovdy,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],rho_nox,v,deltay,0.f);
    first_deriv_product(drhowdz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],rho_nox,w,deltaz,0.f);
    drhodt = -drhoudx - drhovdy - drhowdz;//drhodt = -(node->rho_nox*(dudx+dvdy+dwdz)+drhoudx*node->u+drhovdy*node->v+drhowdz*node->w);
    next_node->rho_nox = node->rho_nox + drhodt * deltat;
    if (next_node->rho_nox < 0.f){
        next_node->rho_nox = 0.f;
    }

    // update rho_h2o
    first_deriv_product(drhoudx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],rho_h2o,u,deltax,0.f);
    first_deriv_product(drhovdy,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],rho_h2o,v,deltay,0.f);
    first_deriv_product(drhowdz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],rho_h2o,w,deltaz,0.f);
    drhodt = -drhoudx - drhovdy - drhowdz;//drhodt = -(node->rho_h2o*(dudx+dvdy+dwdz)+drhoudx*node->u+drhovdy*node->v+drhowdz*node->w);
    next_node->rho_h2o = node->rho_h2o + drhodt * deltat;
    if (next_node->rho_h2o < 0.f){
        next_node->rho_h2o = 0.f;
    }

    float new_rho = next_node->rho_o2+next_node->rho_n2+next_node->rho_fuel+next_node->rho_co2+next_node->rho_nox+next_node->rho_h2o;

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

    float dTxxdx = (2.f*d2udx2 - d2vdxdy - d2wdxdz) * (2.f / 3.f) * node->viscosity;
    float dTyydy = (2.f*d2vdy2 - d2udxdy - d2wdydz) * (2.f / 3.f) * node->viscosity;
    float dTzzdz = (2.f*d2wdz2 - d2udxdz - d2vdydz) * (2.f / 3.f) * node->viscosity;
    float dTxydx = (d2vdx2 + d2udxdy) * node->viscosity;
    float dTxydy = (d2udy2 + d2vdxdy) * node->viscosity;
    float dTyzdy = (d2wdy2 + d2vdydz) * node->viscosity;
    float dTyzdz = (d2vdz2 + d2wdydz) * node->viscosity;
    float dTxzdx = (d2wdx2 + d2udxdz) * node->viscosity;
    float dTxzdz = (d2udz2 + d2wdxdz) * node->viscosity;

    float dpdx, dpdy, dpdz;
    first_deriv(dpdx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],pressure,deltax,node->pressure);
    first_deriv(dpdy,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],pressure,deltay,node->pressure);
    first_deriv(dpdz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],pressure,deltaz,node->pressure);

    float dudt = (-dpdx + dTxxdx + dTxydy + dTxzdz) / rho;
    float dvdt = (-dpdy + dTxydx + dTyydy + dTyzdz) / rho;
    float dwdt = (rho*gravity - dpdz + dTxzdx + dTyzdy + dTzzdz) / rho;

    next_node->u = new_rho > 0.f ? node->u + dudt * deltat : 0.f;
    next_node->v = new_rho > 0.f ? node->v + dvdt * deltat : 0.f;
    next_node->w = new_rho > 0.f ? node->w + dwdt * deltat : 0.f;

    float dV = deltax*deltay*deltaz;
    float dEdx, dEdy, dEdz;
    first_deriv(dEdx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],internal_energy,deltax,boundary_specific_heat*boundary_idling_temp);
    first_deriv(dEdy,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],internal_energy,deltay,boundary_specific_heat*boundary_idling_temp);
    first_deriv(dEdz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],internal_energy,deltaz,boundary_specific_heat*boundary_idling_temp);
    dEdx *= dV;
    dEdy *= dV;
    dEdz *= dV;

    float Txx = (2.f*dudx - dvdy - dwdz) * 2.f/3.f *node->viscosity;
    float Tyy = (2.f*dvdy - dudx - dwdz) * 2.f/3.f *node->viscosity;
    float Tzz = (2.f*dwdz - dudx - dvdy) * 2.f/3.f *node->viscosity;
    float Txy = (dudy + dvdx) * node->viscosity;
    float Txz = (dudz + dwdx) * node->viscosity;
    float Tyz = (dvdz + dwdy) * node->viscosity;
    float dqxdx;
    float dkdx, dTdx, d2Tdx2;
    first_deriv(dkdx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],conductivity,deltax,boundary_conductivity);
    first_deriv(dTdx,nghbrs[nghbrsInd(-1,0,0)],node,nghbrs[nghbrsInd(1,0,0)],temperature,deltax,boundary_idling_temp);
    r0 = nghbrs[nghbrsInd(-1,0,0)];
    r1 = nghbrs[nghbrsInd(0,0,0)];
    r2 = nghbrs[nghbrsInd(1,0,0)];
    d2Tdx2 = (r0 ? r0->temperature : boundary_idling_temp) - 2.f * (r1 ? r1->u : boundary_idling_temp) + (r2 ? r2->u : boundary_idling_temp);
    d2Tdx2 /= (deltax * deltax);
    dqxdx = dkdx * dTdx + node->conductivity * d2Tdx2;

    float dqydy;
    float dkdy, dTdy, d2Tdy2;
    first_deriv(dkdy,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],conductivity,deltay,boundary_conductivity);
    first_deriv(dTdy,nghbrs[nghbrsInd(0,-1,0)],node,nghbrs[nghbrsInd(0,1,0)],temperature,deltay,boundary_idling_temp);
    r0 = nghbrs[nghbrsInd(0,-1,0)];
    r1 = nghbrs[nghbrsInd(0,0,0)];
    r2 = nghbrs[nghbrsInd(0,1,0)];
    d2Tdy2 = (r0 ? r0->temperature : boundary_idling_temp) - 2.f * (r1 ? r1->u : boundary_idling_temp) + (r2 ? r2->u : boundary_idling_temp);
    d2Tdy2 /= (deltay * deltay);
    dqydy = dkdy * dTdy + node->conductivity * d2Tdy2;

    float dqzdz;
    float dkdz, dTdz, d2Tdz2;
    first_deriv(dkdz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],conductivity,deltax,boundary_conductivity);
    first_deriv(dTdz,nghbrs[nghbrsInd(0,0,-1)],node,nghbrs[nghbrsInd(0,0,1)],temperature,deltax,boundary_idling_temp);
    r0 = nghbrs[nghbrsInd(0,0,-1)];
    r1 = nghbrs[nghbrsInd(0,0,0)];
    r2 = nghbrs[nghbrsInd(0,0,1)];
    d2Tdz2 = (r0 ? r0->temperature : boundary_idling_temp) - 2.f * (r1 ? r1->u : boundary_idling_temp) + (r2 ? r2->u : boundary_idling_temp);
    d2udz2 /= (deltaz * deltaz);
    dqzdz = dkdz * dTdz + node->conductivity * d2Tdz2;

    float dEdt = node->dQ / deltat;
    dEdt += rho*gravity*node->w;
    dEdt -= dEdx*node->u + dEdy*node->v + dEdz*node->w + 
            node->internal_energy*dV*(dudx + dvdy + dwdz);
    dEdt -= dpdx*node->u + dpdy*node->v + dpdz*node->w +
            node->pressure * (node->u+node->v+node->w);
    dEdt += dudx*Txx + dudy*Txy + dudz * Txz + 
            node->u * (dTxxdx + dTxydy + dTxzdz);
    dEdt += dvdx*Txy + dvdy*Tyy + dvdz*Tyz +
            node->v * (dTxydx + dTyydy + dTyzdz);
    dEdt += dwdx*Txz + dwdy*Tyz + dwdz*Tzz +
            node->w * (dTxzdx + dTyzdy + dTzzdz);
    dEdt -= dqxdx + dqydy + dqzdz;

    next_node->internal_energy = node->internal_energy + (dEdt / dV) * deltat;
    next_node->temperature = (new_rho <= 0 ? 0.f : 375 +
                              next_node->internal_energy-268000.f/(723.5f*rho));
    float temp = next_node->temperature;
    next_node->viscosity = 0.000001458f * sqrt(temp*temp*temp) / (temp + 110.4f);
    next_node->conductivity = 0.000001458f * sqrt(temp*temp*temp) / (temp + 110.4f);
    //TODO: find a way to get square root in CUDA
    float nV = 0.f;
    nV += next_node->rho_o2 / o2_molar_mass;
    nV += next_node->rho_n2 / n2_molar_mass;
    nV += next_node->rho_fuel / fuel_molar_mass;
    nV += next_node->rho_co2 / co2_molar_mass;
    nV += next_node->rho_nox / nox_molar_mass;
    nV += next_node->rho_h2o / h2o_molar_mass;
    next_node->pressure = (nV <= 0.f ? 0.f : nV * 8.31 * temp);
}

__global__ void kernelSpark(){
    int nodeX = blockIdx.x * blockDim.x + threadIdx.x;
    int nodeY = blockIdx.y * blockDim.y + threadIdx.y;
    int nodeZ = blockIdx.z * blockDim.z + threadIdx.z;
    if (nodeX >= cuImageData.nodeLength || nodeY >= cuImageData.nodeWidth
                    || nodeZ >= cuImageData.nodeDepth){
        return;
    }
    int nodeIdx = getNodeIdx(nodeX, nodeY, nodeZ);
    if (!((cuImageData.nodes[nodeIdx].temperature >= fuel_autoignition_point) || 
            (cuImageData.nodes[nodeIdx].temperature >= fuel_flash_point &&
            false))){
            // sparks.find(index) != sparks.end()))){ 
                //TODO: check if ^this works
        return;
    }
    float rho_o2 = cuImageData.nodes[nodeIdx].rho_o2;
    float nV_o2 = rho_o2 / o2_molar_mass;
    float rho_fuel = cuImageData.nodes[nodeIdx].rho_fuel;
    float nV_fuel = rho_fuel / fuel_molar_mass;

    float delta_o2, delta_n2, delta_fuel, delta_co2, delta_nox, delta_h2o;
    if (2.f*nV_o2 >= 25.f*nV_fuel){
        // reaction is limited by fuel
        delta_fuel = -reaction_rate_coefficient * nV_fuel * cuImageData.dt;
        delta_o2 = -reaction_rate_coefficient * nV_o2 * cuImageData.dt;
        delta_nox = -(2.f*delta_o2 - 25.f*delta_fuel);
    }
    else{
        // reaction is limited by air
        delta_o2 = -reaction_rate_coefficient * nV_o2 * cuImageData.dt;
        delta_fuel = delta_o2 / 12.5f;
        delta_nox = 0.f;
    }
    delta_co2 = -delta_fuel * 8.f;
    delta_h2o = -delta_fuel * 9.f;
    delta_n2 = -delta_nox / 2.f;

    cuImageData.nodes[nodeIdx].rho_o2 += delta_o2 * o2_molar_mass;
    cuImageData.nodes[nodeIdx].rho_n2 += delta_n2 * n2_molar_mass;
    cuImageData.nodes[nodeIdx].rho_fuel += delta_fuel * fuel_molar_mass;
    cuImageData.nodes[nodeIdx].rho_co2 += delta_co2 * co2_molar_mass;
    cuImageData.nodes[nodeIdx].rho_nox += delta_nox * nox_molar_mass;
    cuImageData.nodes[nodeIdx].rho_h2o += delta_h2o * h2o_molar_mass;

    cuImageData.nodes[nodeIdx].dQ = -delta_o2 * o2_formation_enthalpy
                        -delta_n2 * n2_formation_enthalpy
                        -delta_fuel * fuel_formation_enthalpy
                        -delta_co2 * co2_formation_enthalpy
                        -delta_nox * nox_formation_enthalpy
                        -delta_h2o * h2o_formation_enthalpy;
}

__device__ __inline__ float4 shadePixel(int nodeIdx, float4 pixel){
    Node curr_node = cuImageData.nodes[nodeIdx];
    float pressure = curr_node.pressure;
    float temperature = curr_node.temperature;
    float internal_energy = curr_node.internal_energy;
    float r = (pressure < PRESSURE_LOWER_BOUND ? 0.f : 
               ((pressure > PRESSURE_UPPER_BOUND) ? 255.f :
               (255.f * (pressure - PRESSURE_LOWER_BOUND))));
    float g = (temperature < TEMP_LOWER_BOUND ? 0 : 
               (temperature > TEMP_UPPER_BOUND ? 255.f :
               (255.f * (temperature - TEMP_LOWER_BOUND))));
    float b = (internal_energy < ENERGY_LOWER_BOUND ? 0 : 
               (internal_energy > ENERGY_UPPER_BOUND ? 255.f :
               (255.f * (internal_energy - ENERGY_LOWER_BOUND))));
    //TODO: find a way to handle multiple layers.
    return make_float4(r, g, b, 255.f);
}

__global__ void kernelCopyNodes(){
    int nodeX = blockIdx.x * blockDim.x + threadIdx.x;
    int nodeY = blockIdx.y * blockDim.y + threadIdx.y;
    int nodeZ = blockIdx.z * blockDim.z + threadIdx.z;
    if (nodeX >= cuImageData.nodeLength || nodeY >= cuImageData.nodeWidth
                    || nodeZ >= cuImageData.nodeDepth){
        return;
    }
    int nodeIdx = getNodeIdx(nodeX, nodeY, nodeZ);
    cuImageData.nodes[nodeIdx] = cuImageData.newNodes[nodeIdx];
}

__global__ void kernelRenderImage(){
    int imageX = blockIdx.x * blockDim.x + threadIdx.x;
    int imageY = blockIdx.y * blockDim.y + threadIdx.y;

    int width = cuImageData.imageWidth;
    int height = cuImageData.imageHeight;

    if (imageX >= width || imageY >= height) return;

    float4* imgPtr = (float4*) (&cuImageData.imageData[4 * (imageY * width + imageX)]);

    float4 pixel = *imgPtr;

    int nodeIdx;
    for (int z = 0; z < cuImageData.nodeDepth; z++){
        nodeIdx = imageX + imageY * cuImageData.nodeLength +
                    z * cuImageData.nodeLength * cuImageData.nodeWidth;
        pixel = shadePixel(nodeIdx, pixel);
    }
    *imgPtr = pixel;

    return;
}

// from hw2 assignment
__global__ void kernelClearImage(float r, float g, float b, float a){
    int imageX = blockIdx.x * blockDim.x + threadIdx.x;
    int imageY = blockIdx.y * blockDim.y + threadIdx.y;

    int width = cuImageData.imageWidth;
    int height = cuImageData.imageHeight;

    if (imageX >= width || imageY >= height)
        return;

    int offset = 4 * (imageY * width + imageX);
    float4 value = make_float4(r, g, b, a);

    // write to global memory: As an optimization, I use a float4
    // store, that results in more efficient code than if I coded this
    // up as four seperate fp32 stores.
    *(float4*)(&cuImageData.imageData[offset]) = value;
}

void CudaVisualizer::simulateSteps(){
    /*
     * uh i don't think nodes are going to fit into the warps here. Like I'm
     * pretty sure we'll have to move the properties to some array and do
     * with that instead of just passing the structs into the device memory.
     * It might just be easier to do just the visualizer on this. 
     */
    dim3 blockDim(16, 16, 4);
    dim3 gridDim((nodeWidth + 15)/16,
                 (nodeLength + 15)/16,
                 (nodeDepth + 3)/4);
    for(int i = 0; i < numIterations; i++){
        double startTime = CycleTimer::currentSeconds();
        kernelSpark<<<blockDim, gridDim>>>();
        kernelSimStep<<<blockDim, gridDim>>>();
        double endTime = CycleTimer::currentSeconds();
        kernelCopyNodes<<<blockDim, gridDim>>>();
        cudaError_t errCode = cudaPeekAtLastError();
        if (errCode != cudaSuccess) {
            fprintf(stderr, "WARNING: A CUDA error occured on iteration %d: code=%d, %s\n", i, errCode, cudaGetErrorString(errCode));
        }
        printf("iteration %d took %f seconds on CUDA\n", i, endTime-startTime);
    }
    printf("%d\n", nodeSize);
    cudaMemcpy(nodes, cuNewNodes, nodeSize * sizeof(Node),
                cudaMemcpyDeviceToHost);
    return;
}

std::vector<Node>
CudaVisualizer::getNodes(){
    const std::vector<Node> res(nodes, nodes+nodeSize);
    return res;
}

void
CudaVisualizer::allocOutputImage(){
    if (image) delete image;
    image = new Image(nodeWidth, nodeLength);
    return;
}

void CudaVisualizer::clearImage(){
    dim3 blockDim(16, 16, 1);
    dim3 gridDim(
        (image->width + blockDim.x - 1) / blockDim.x,
        (image->height + blockDim.y - 1) / blockDim.y);
    kernelClearImage<<<gridDim, blockDim>>>(1.f, 1.f, 1.f, 1.f);
    return;
}

void
CudaVisualizer::shade(){
    int imageWidth = image->width;
    int imageHeight = image->height;
    dim3 blockDim(32, 32);
    dim3 chunkDim((imageWidth + 31) / 32,
                   (imageHeight + 31) / 32);
    kernelRenderImage<<<chunkDim, blockDim>>>();
    return;
}

CudaVisualizer::CudaVisualizer(){
    image = NULL;
    nodes = NULL;
    sparks = NULL;
    
    sparkLen = 0;
    nodeWidth = 0;
    nodeLength = 0;
    nodeDepth = 0;
    nodeSize = 0;
    dy = 0;
    dz = 0;
    dt = 0;
    dx = 0;

    cuImage = NULL;
    cuNodes = NULL;
    cuNewNodes = NULL;
}

CudaVisualizer::~CudaVisualizer(){
    if (image) delete image;
    if (nodes) delete nodes;
    if (cuImage){
        cudaFree(cuImage);
        cudaFree(cuNodes);
        cudaFree(cuNewNodes);
    }
}

void
CudaVisualizer::init(){
    if (cuImage){
        cudaFree(cuImage);
        cudaFree(cuNodes);
        cudaFree(cuNewNodes);
        cudaFree(cuSparks);
    }
    if (!image) allocOutputImage();
    cudaMalloc(&cuImage, 4 * sizeof(float) * image->width * image->height);
    cudaMalloc(&cuNodes, sizeof(Node) * nodeSize);
    cudaMalloc(&cuNewNodes, sizeof(Node) * nodeSize);
    cudaMalloc(&cuSparks, sizeof(int) * sparkLen);

    GlobalConstants params;
    params.imageWidth = image->width;
    params.imageHeight = image->height;
    params.imageData = cuImage;
    params.nodeWidth = nodeWidth;
    params.nodeLength = nodeLength;
    params.nodeDepth = nodeDepth;
    params.sparkLen = sparkLen;
    params.nodes = cuNodes;
    params.newNodes = cuNewNodes;
    params.dx = dx;
    params.dy = dy;
    params.dz = dz;
    params.dt = dt;

    cudaMemcpyToSymbol(cuImageData, &params, sizeof(GlobalConstants));
    cudaMemcpy(cuNodes, nodes, sizeof(Node) * nodeSize,
                cudaMemcpyHostToDevice);
}

const Image*
CudaVisualizer::getImage(){
    cudaMemcpy(
        image->data, cuImage, sizeof(float) * 4 * image->width * image->height,
        cudaMemcpyDeviceToHost
    );
    return image;
}

void
CudaVisualizer::setParams(std::vector<Node>& new_nodes,
                          std::vector<int> spark_vec,
                          const stepParams params,
                          int iterations, int sparkSize){
    nodes = &new_nodes[0];
    nodeSize = new_nodes.size();
    nodeWidth = params.width;
    nodeLength = params.length;
    nodeDepth = params.depth;
    dx = params.deltax;
    dy = params.deltay;
    dz = params.deltaz;
    dt = params.deltat;
    sparkLen = sparkSize;
    sparks = &spark_vec[0];
    numIterations = iterations;
}


void 
CudaVisualizer::render(){
    if (!image){
        allocOutputImage();
    }
    init();
    clearImage();
    double startTime = CycleTimer::currentSeconds();
    shade();
    double endTime = CycleTimer::currentSeconds();
    printf("shading took %f seconds\n", endTime-startTime);
    writePPMImage(getImage(), "imageOutput.ppm");
    return;
}
