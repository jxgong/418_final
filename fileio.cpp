#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "common.h"

std::vector<Node> loadFromFile(std::string filename, stepParams *params, int *num_iterations){
    std::ifstream inFile;
    inFile.open(filename);
    std::string line;
    std::vector<Node> result;
    if (!inFile) return result;
//    TODO: randomly generate node info if no file found
    inFile >> params->deltax;
    inFile >> params->deltay;
    inFile >> params->deltaz;
    inFile >> params->deltat;
    inFile >> *num_iterations;
    inFile >> params->length;
    inFile >> params->width;
    inFile >> params->depth;

    while (std::getline(inFile, line))
    {
        Node node;
        std::stringstream sstream(line);
        std::string str;
        std::getline(sstream, str, ' ');
        node.rho_o2 = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.rho_n2 = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.rho_fuel = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.rho_co2 = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.rho_nox = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.rho_h2o = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.pressure = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.temperature = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.viscosity = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.internal_energy = (float)atof(str.c_str());
        std::getline(sstream, str, ' ');
        node.conductivity = (float)atof(str.c_str());
        result.push_back(node);
    }
    return result;
}

bool saveToFile(std::vector<Node> data, std::string filename){
    std::ofstream file(filename);
    if (!file)
    {
        std::cout << "error writing file \"" << filename << "\"" << std::endl;
        return false;
    }
    file << std::setprecision(9);
    Node curr_node;
    for (unsigned int i = 0; i < data.size(); i++){
        curr_node = data.at(i);
        file << curr_node.rho_o2 << " " << curr_node.rho_n2 << " "
            << curr_node.rho_fuel << " " << curr_node.rho_co2 << " " 
            << curr_node.rho_nox << " " << curr_node.rho_h2o << " " 
            << curr_node.pressure << " " << curr_node.temperature << " " 
            << curr_node.viscosity << " " << curr_node.internal_energy << " "
            << curr_node.conductivity << std::endl;
    }
    file.close();
    return true;
}
