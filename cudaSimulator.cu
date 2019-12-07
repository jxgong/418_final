#include "common.h"
#include "physics.h"
#include <vector>
#include "CycleTimer.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>

__global__ void kernelSimSteps(){
    return;
}

void simulateStepCuda(std::vector<Node>& new_nodes,
                      std::vector<Node>& nodes,
                      const stepParams params){
    
    return;
}