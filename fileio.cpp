#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "common.h"

std::vector<Node> loadFromFile(std::string filename){
    std::ifstream inFile;
    inFile.open(filename);
    std::string line;
    std::vector<Node> result;
    if (!inFile) return result;
//    TODO: randomly generate node info if no file found
    while (std::getline(inFile, line))
    {
        Node node;
        std::stringstream sstream(line);
        std::string str;
        std::getline(sstream, str, ' ');
        node.rho_air = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.rho_fuel = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.rho_co2 = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.rho_nox = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.pressure = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.temperature = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.viscosity = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.u = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.v = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.w = (float)atof(str.c_str());
        result.push_back(node);
    }
    return result;
}

bool saveToFile(std::vector<Node> data, std::string filename){
    return true;
}