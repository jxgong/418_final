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
}

__global__ void kernelSimSteps(){
    return;
}
 __global__ void kernelRenderImage(){
    return;
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
CudaVisualizer::render(std::vector<Node>& nodes,
                        const stepParams params){
    if (!image){
        allocOutputImage(params.length, params.width);
    }
    clearImage();
    return;
}
