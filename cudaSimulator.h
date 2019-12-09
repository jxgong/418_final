#ifndef CUDA_SIMULATOR_H
#define CUDA_SIMULATOR_H

#include "image.h"
#include "common.h"

class CudaVisualizer {
    private:
        Image* image;
        Node* nodes;
        int* sparks;
        int sparkLen;
        int nodeWidth;
        int nodeLength;
        int nodeDepth;
        float dx, dy, dz, dt;
        int nodeSize;
        int numIterations;

        float* cuImage;
        Node* cuNodes;
        Node* cuNewNodes;
        int* cuSparks;

    public:
        CudaVisualizer();
        virtual ~CudaVisualizer();
        const Image* getImage();
        void allocOutputImage();
        void clearImage();
        void simulateSteps();
        void render();
        void setParams(std::vector<Node>& nodes, std::vector<int> spark_vec,
                        const stepParams params,
                        int iterations, int sparkSize);
        void shade();
        void init();
        std::vector<Node> getNodes();
};

#endif