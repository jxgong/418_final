#include <string>

class Node
{
    public:
        float m_air, m_fuel, m_co2, m_nox, temp, u, v, w;
};


struct StartupOptions
{
    int numIterations = 1;
    float length, width, height;
    std::string inputFile = "input.txt";
    bool checkCorrectness = true;
};

StartupOptions parseOptions(int argc, char *argv[]);