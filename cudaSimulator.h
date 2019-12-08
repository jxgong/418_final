#ifndef CUDA_SIMULATOR_H
#define CUDA_SIMULATOR_H

#include "image.h"

class CudaVisualizer {
    private:
        Image* image;

    public:
        const Image* getImage();
        void allocOutputImage();
        void clearImage();
        void render();
};

#endif