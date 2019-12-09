#include "common.h"
#include "physics.h"
#include <vector>
#include "image.h"
#include "CycleTimer.h"
#include "cudaSimulator.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>


__constant__ int PRESSURE_LOWER_BOUND = 0.f;
__constant__ int PRESSURE_UPPER_BOUND = 100.f;
__constant__ int TEMP_UPPER_BOUND = 100.f;
__constant__ int TEMP_LOWER_BOUND = 0.f;
__constant__ int ENERGY_LOWER_BOUND = 0.f;
__constant__ int ENERGY_UPPER_BOUND = 100.f;
//TODO: check if these bounds are reasonable.

struct GlobalConstants{
    int imageWidth;
    int imageHeight;
    float* imageData;
    
    int nodeWidth, nodeLength, nodeDepth;
    float* pressure;
    float* temperature;
    float* internal_energy;
};

__constant__ GlobalConstants cuImageData;

__global__ void kernelSimSteps(){
    return;
}

__device__ __inline__ float4 shadePixel(int nodeIdx, float4 pixel){
    float pressure = cuImageData.pressure[nodeIdx];
    float temperature = cuImageData.temperature[nodeIdx];
    float internal_energy = cuImageData.internal_energy[nodeIdx];
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
        //TODO: check if i'm mixing up length and width here.
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

void simulateStepCuda(std::vector<Node>& new_nodes,
                      std::vector<Node>& nodes,
                      const stepParams params){
    /*
     * uh i don't think nodes are going to fit into the warps here. Like I'm
     * pretty sure we'll have to move the properties to some array and do
     * with that instead of just passing the structs into the device memory.
     * It might just be easier to do just the visualizer on this. 
     */
    return;
}
void
CudaVisualizer::allocOutputImage(int length, int width){
    if (image) delete image;
    image = new Image(width, length);
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
    dim3 chunkDim((imageWidth + 32) / 32,
                   (imageHeight + 32) / 32);
    kernelRenderImage<<<chunkDim, blockDim>>>();
    return;
}


void
CudaVisualizer::init(){
    cudaMalloc(&cuImage, 4 * sizeof(float) * image->width * image->height);
    cudaMalloc(&cuPressure, sizeof(float) * nodes.size());
    cudaMalloc(&cuTemp, sizeof(float) * nodes.size());
    cudaMalloc(&cuInternalEnergy, sizeof(float) * nodes.size());

    GlobalConstants params;
    params.imageWidth = image->width;
    params.imageHeight = image->height;
    params.imageData = cuImage;
    params.pressure = cuPressure;
    params.temperature = cuTemp;
    params.internal_energy = cuInternalEnergy;
    params.nodeWidth = step_params.width;
    params.nodeLength = step_params.length;
    params.nodeDepth = step_params.depth;

    cudaMemcpyToSymbol(cuImageData, &params, sizeof(GlobalConstants));
}
void 
CudaVisualizer::render(std::vector<Node>& nodes,
                        const stepParams params){
    if (!image){
        allocOutputImage(params.length, params.width);
    }
    clearImage();
    shade();
    return;
}
