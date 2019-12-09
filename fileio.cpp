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

    for (int i = 0; i < params->length * params->width * params->depth; i++){
        Node node;
        inFile >> node.rho_o2;
        inFile >> node.rho_n2;
        inFile >> node.rho_fuel;
        inFile >> node.rho_co2;
        inFile >> node.rho_nox;
        inFile >> node.rho_h2o;
        inFile >> node.pressure;
        inFile >> node.temperature;
        inFile >> node.viscosity;
        inFile >> node.internal_energy;
        inFile >> node.conductivity;
        inFile >> node.u;
        inFile >> node.v;
        inFile >> node.w;
        node.dQdt = 0.f;
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
