#ifndef CUDA_SIMULATOR_H
#define CUDA_SIMULATOR_H

#include "image.h"
#include "common.h"

class CudaVisualizer {
    private:
        Image* image;
        Node* nodes;
        int nodeWidth;
        int nodeLength;
        int nodeDepth;
        int nodeSize;

        float* cuImage;
        float* cuPressure;
        float* cuTemp;
        float* cuInternalEnergy;
        Node* cuNodes;

    public:
        CudaVisualizer();
        ~CudaVisualizer();
        const Image* getImage();
        void allocOutputImage();
        void clearImage();
        void render(std::vector<Node>& nodes, const stepParams params);
        void setParams(std::vector<Node>& nodes, const stepParams params);
        void shade();
        void init();
};

#endif