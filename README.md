# Title
#### Jason Gong, Fernando Melean
## Project URL: https://jxgong.github.io/fuel_ignition_simulator/
## Summary:
We are going to parallelize simulations of fuel ignition and flame propagation simulate using a GPU.
## Background:
The field of fluid dynamics involves a large amount of computation. In order to properly simulate something like fuel ignition, pressure, temperature, and the velocity of the fluid all have an effect. Although a fluid is a continuous substance, computers cannot simulate continuous things, so we need to find a way to approximate a real world fluid in a discrete system. Instead of having a fluid as a continuous object, we can simulate each particle in parallel.
Fuel ignition and propagation is a chemical interaction dependent on the above factors. This reaction has effects on the surrounding particles as energy gets released, which will affect further iterations of the simulation.
## The Challenge:
There are many different parts of this problem that can be parallelized. As mentioned above, each particle is affected by pressure and temperature, which we'll have to keep track of. Also, since ignition is a chemical reaction, we'll have to keep track of what's available to simulate it as well. The equations to discretize the fluids are also decently complex. Also, calculating the states of the particles depend on each other, so we'll have to find a way to find dependencies between particles.
## Goals and Deliverables:
### Platform Choice:
We'll use the GHC cluster machines to run our simulations. The GHC cluster computers have 8 Intel processors and NVIDIA GeForce GTX 1080 GPUs
## Schedule:
| Week  | Task                                         |
|-------|----------------------------------------------|
| 10/27 | Finish Project Proposal                      |
| 11/3  | Sequential Simulator with Fuel ignition      |
| 11/10 | Naive Parallel Sim for updating particles    |
| 11/17 | Optimize Parallelism, start doing visualizer |
| 11/24 | Finish visualizer                            |
| 12/1  | Polish                                       |
| 12/8  | Poster                                       |
