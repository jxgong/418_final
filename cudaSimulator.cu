#include "common.h"
#include "physics.h"
#include <vector>
#include "image.h"
#include "CycleTimer.h"
#include "cudaSimulator.h"
#include "ppm.h"

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
    
    int nodeWidth, nodeLength, nodeDepth;
    float dx, dy, dz, dt;
    Node* nodes;
    Node* newNodes;
};

__constant__ GlobalConstants cuImageData;

__global__ void kernelSimStep(){
    int nodeX = blockIdx.x * blockDim.x + threadIdx.x;
    int nodeY = blockIdx.y * blockDim.y + threadIdx.y;
    int nodeZ = blockIdx.z * blockDim.z + threadIdx.z;
    int nodeIdx = nodeX + params.nodeLength * (nodeY + nodeWidth * nodeZ);
    if (!((params.nodes[index].temperature >= fuel_autoignition_point) || 
            (params.nodes[index].temperature >= fuel_flash_point &&
            sparks.find(index) != sparks.end()))){ 
                //TODO: check if ^this works
        return;
    }
    float rho_o2 = params.nodes[index].rho_o2;
    float nV_o2 = rho_o2 / o2_molar_mass;
    float rho_fuel = params.nodes[index].rho_fuel;
    float nV_fuel = rho_fuel / fuel_molar_mass;

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

    params.nodes[index].rho_o2 += delta_o2 * o2_molar_mass;
    params.nodes[index].rho_n2 += delta_n2 * n2_molar_mass;
    params.nodes[index].rho_fuel += delta_fuel * fuel_molar_mass;
    params.nodes[index].rho_co2 += delta_co2 * co2_molar_mass;
    params.nodes[index].rho_nox += delta_nox * nox_molar_mass;
    params.nodes[index].rho_h2o += delta_h2o * h2o_molar_mass;

    params.nodes[index].dQdt = -delta_o2 * o2_formation_enthalpy
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
        kernelSimStep<<<blockDim, gridDim>>>();
        cudaMemcpy(cuNodes, cuNewNodes,
                   nodeSize * sizeof(Node),
                   cudaMemcpyDeviceToDevice);
        double endTime = CycleTimer::currentSeconds();
        printf("iteration %d took %f seconds on CUDA\n", i, endTime-startTime);
    }
    return;
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
    }
    cudaMalloc(&cuImage, 4 * sizeof(float) * image->width * image->height);
    cudaMalloc(&cuNodes, sizeof(Node) * nodeSize);
    cudaMalloc(&cuNewNodes, sizeof(Node) * nodeSize);

    GlobalConstants params;
    params.imageWidth = image->width;
    params.imageHeight = image->height;
    params.imageData = cuImage;
    params.nodeWidth = nodeWidth;
    params.nodeLength = nodeLength;
    params.nodeDepth = nodeDepth;
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
                          const stepParams params,
                          int iterations){
    nodes = &new_nodes[0];
    nodeSize = new_nodes.size();
    nodeWidth = params.width;
    nodeLength = params.length;
    nodeDepth = params.depth;
    dx = params.deltax;
    dy = params.deltay;
    dz = params.deltaz;
    dt = params.deltat;
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
