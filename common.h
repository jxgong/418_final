#include <string>
#include <vector>

class Node
{ 
    public:
        // physical properties
        float rho_air, rho_fuel, rho_co2, rho_nox;
        float pressure, temperature;
        float viscosity, internal_energy; 
        // dynamic properties
        float u, v, w; 
        float get_rho();

};


struct StartupOptions
{
    int numIterations = 1;
    float length, width, height;
    std::string inputFile = "input.txt";
    bool checkCorrectness = true;
};

class stepParams{
    public:
        int length, width, depth;
        float deltat, deltax, deltay, deltaz;
        std::vector<int> sparks;
        std::vector<Node> valves;
};

StartupOptions parseOptions(int argc, char *argv[]);
