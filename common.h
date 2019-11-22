#include <string>

class Node
{
    public:
        float rho_air, rho_fuel, rho_co2, rho_nox, pressure, temperature, viscosity; // physical properties
        float u, v, w; // dynamic properties
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