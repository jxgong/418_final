# Fuel Ignition And Flame Propagation CFD
#### Jason Gong, Fernando Melean
## Project URL: 
https://jxgong.github.io/fuel_ignition_simulator/
## Summary:
We are going to parallelize simulations of fuel ignition and flame propagation simulate using a GPU.
## Background:
CFD (computational fluid-dynamic) models are used to study fluid motion for a variety of applications, ranging from medicine to propulsion to water filtering. Simulation of fuel ignition and flame propagation with CFDs are of particular interest to those working with propulsion systems such as cars, rockets, or planes in studying how to improve fuel efficiency and reduce pollutant production in combustion engines. However, these models are very computationally expensive, often taking days for a single simulation, limiting exploration of design space. Our goal will be to speed-up common operations for these models such as advection and combustion using GPUs. 

These calculations will rely on applications of Navier-Stokes differential equations and heat transfer equations. Since we cannot simulate this analytically, we will instead use the finite element method, breaking down our volume of interest into many small finite volume elements, and apply finite difference methods (FDMs) such as taylor series to approximate physical and chemical properties of each element including temperature, density, fuel-to-air ratio, etc. Using these FDMs will have a few effects, (1) distant regions in the volume will have negligible effects on each other allowing them to be performed in parallel and (2) operations on the elements will consist of a few, highly repetitive operations making them good candidates for SIMD speedup.
## The Challenge:
There are many different parts of this problem that can be parallelized. As mentioned above, each finite element is affected by physical properties such as temperature and density, which we'll have to keep track of. Also, since ignition is a chemical reaction, we'll have to keep track of each element’s chemical composition (e.g. air-to-fuel ratio) as well to simulate them. The equations to discretize the fluids are also decently complex, relying heavily on taylor series and other finite difference methods, which introduce dependencies between the finite elements as each element will be updated based on many nearby elements. However, the dependencies reach farther than immediate neighbors due to higher complexity of the differential equation. This, however, suggests it will have good locality for memory accesses. The main divergence we anticipate in the workload is that at only a fraction of cells may be burning fuel leading to those cells requiring extra computation relative to cells only experiencing fluid motion.  
## Goals and Deliverables:
We plan to achieve a speedup over the sequential code on CPU on the order of 1000’s for the parallel program on GPU. We believe we can do this because (1) the computation is very intensive, yet very similar for each element, making it a good candidate for SPMD, (2) has good locality for memory accesses since elements only reference nearby elements so the computation is less likely to be memory bound, and (3) the cluster GPUs have a theoretical peak speedup of 81,290x due to having 2560 cores with 32-wide warps. 

If ahead of schedule, we hope to achieve high speedup for very large simulations that might not fit on a single machine. This might involve using MPI to divide regions of the simulation across machines and send boundary elements between them. Since computation relies heavily on nearby elements, the communication to computation ratio for this should be fairly low. 

Our demo should include our visualizer of the fuel ignition. We will also show speedups graphs between the sequential and parallel versions of the code for each of the important operations such as advection and combustion.
### Deliverables:
* Sequential Implementation of CFD for fuel ignition and flame propagation
* Parallel Implementation of CFD on cluster GPUs using CUDA
* Speedup graphs for advection and combustion
* Visualization of ignition and flame propagation
### Platform Choice:
We'll use the GHC cluster machines to run our simulations. The GHC cluster computers have 8 Intel processors and NVIDIA GeForce GTX 1080 GPUs
## Schedule:
![alt text](https://raw.githubusercontent.com/jxgong/fuel_ignition_simulator/master/schedulev1.jpg)

## Checkpoint 2019-11-21
Progress on the project has been disappointing so far. Fluid dynamics is a complicated topic and learning the equations needed takes a lot of mathematical and chemistry background that was difficult to understand and took a long time to learn. The progress we did make was mostly sequential and we haven't done much optimization.

We've managed to implement basic conservation of mass and conservation of momentum in our simulate step. Most of our other progress was based around understanding fluid dynamics. Despite this setback, our deliverables should remain largely the same. Specifically, we will still get a sequential and parallel implementation of the CFD, and we'll still present a few speedup graphs for the speedup. However, we may not have time for the visualizer, and we may not be able to achieve as much speedup as we initially desired.

We plan to show speedup graphs at the speedup session, and if we manage to get a basic visualizer working, we may be able to show a demo as well.

Concerns at this point include time management, the fact that I (Jason) still haven't really gotten the concepts needed for fluid dynamics.

## Update 2019-12-09

[Final project is here](https://github.com/jxgong/fuel_ignition_simulator/raw/master/Fuel%20Ignition%20Sim%20Report.pdf)
