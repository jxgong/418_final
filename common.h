
#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <vector>

#define first_deriv(varname,prev,curr,next,field,delta, default_value)\
                if (prev && next){ \
                    varname = ((next)->field - (prev)->field) / (2.f*(delta)); \
                } \
                else if (prev){ \
                    varname = (2.f*(default_value) - (curr)->field - (prev)->field) / (2.f*(delta));\
                } \
                else if (next){ \
                    varname = ((next)->field + (curr)->field - 2.f*(default_value)) / (2.f*(delta));\
                } \
                else{ \
                    varname = 0.f; \
                }

#define dim2Index(x,y,z,Lx,Ly) ((x) + ((y)*(Lx)) + ((z)*(Lx)*(Ly)))
#define nghbrsInd(x,y,z) (((x)+1) + 3 * ((y)+1) + 9 * ((z)+1))

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
        // std::unordered_set<int> sparks;
};

StartupOptions parseOptions(int argc, char *argv[]);

#endif