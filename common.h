
#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <vector>

class Node
{ 
    public:
        // physical properties
        float rho_o2, rho_n2, rho_fuel, rho_co2, rho_nox, rho_h2o;
        float pressure, temperature;
        float viscosity, internal_energy; 
        float conductivity;
        // dynamic properties
        float u, v, w; 
        // thermal properties
        float dQdt;
        float get_rho();

};


struct StartupOptions
{
    int numIterations;
    float length, width, height;
    std::string inputFile;
    bool checkCorrectness;
};

class stepParams{
    public:
        int length, width, depth;
        float deltat, deltax, deltay, deltaz;
        std::vector<int> sparks;
        std::vector<Node> valves;
};

StartupOptions parseOptions(int argc, char *argv[]);

#endif
