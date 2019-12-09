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
        float dx, dy, dz, dt;
        int nodeSize;
        int numIterations;

        float* cuImage;
        Node* cuNodes;
        Node* cuNewNodes;

    public:
        CudaVisualizer();
        virtual ~CudaVisualizer();
        const Image* getImage();
        void allocOutputImage();
        void clearImage();
        void simulateSteps();
        void render();
        void setParams(std::vector<Node>& nodes, const stepParams params,
                        int iterations);
        void shade();
        void init();
};

#endif