#include "common.h"
#include "physics.h"
#include <vector>
#include "image.h"
#include "CycleTimer.h"
#include "cudaSimulator.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>

struct GlobalConstants{
    int imageWidth;
    int imageHeight;
    float* imageData;

    float* pressure;
    float* temperature;
    float* internal_energy;
};

__constant__ GlobalConstants cuImageData;

__global__ void kernelSimSteps(){
    return;
}

__global__ void kernelRenderImage(){
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
    return;
}

void
CudaVisualizer::shade(){
    
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
