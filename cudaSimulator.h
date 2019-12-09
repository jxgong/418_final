#ifndef CUDA_SIMULATOR_H
#define CUDA_SIMULATOR_H

#include "image.h"
#include "common.h"

class CudaVisualizer {
    private:
        Image* image;
        std::vector<Node>& nodes;
        const stepParams step_params;

        float* cuImage;
        float* cuPressure;
        float* cuTemp;
        float* cuInternalEnergy;

    public:
        const Image* getImage();
        void allocOutputImage(int width, int height);
        void clearImage();
        void render(std::vector<Node>& nodes, const stepParams params);
        void shade();
        void init();
};

#endif