#include <iostream>
#include <fstream>
#include "physics.h"
#include <cstdlib>

const char *filename = "input/uniform_100-100-100.txt";
float deltax = 0.001f;
float deltay = 0.001f;
float deltaz = 0.001f;
float deltat = 0.05f;
int numIterations = 600;
int length = 100;
int width = 100;
int depth = 100;

static float rho_o2(int x, int y, int z){
    return 0.2854f;
}
static float rho_n2(int x, int y, int z){
    return 0.9396f;
}
static float rho_fuel(int x, int y, int z){
    return 1.225f / 13.6f;
}
static float rho_co2(int x, int y, int z){
    return 0.f;
}
static float rho_nox(int x, int y, int z){
    return 0.f;
}
static float rho_h2o(int x, int y, int z){
    return 0.f;
}
static float temperature(int x, int y, int z){
    return boundary_idling_temp;
}
static float u(int x, int y, int z){
    return 0.2f * (((float) rand()) / (float) RAND_MAX) - 0.1f;
}
static float v(int x, int y, int z){
    return 0.2f * (((float) rand()) / (float) RAND_MAX) - 0.1f;
}
static float w(int x, int y, int z){
    return 0.2f * (((float) rand()) / (float) RAND_MAX) - 0.1f;
}

int main(int argc, char const *argv[])
{
    std::ofstream file;
    file.open(filename);

    file << deltax << std::endl <<
            deltay << std::endl <<
            deltaz << std::endl <<
            deltat << std::endl;
    file << numIterations << std::endl;
    file << length << std::endl <<
            width << std::endl <<
            depth << std::endl;

    for (int k = 0; k < depth; k++){
        for (int j = 0; j < width; j++){
            for (int i = 0; i < length; i++){
                Node node;
                node.rho_o2 = rho_o2(i,j,k);
                node.rho_n2 = rho_n2(i,j,k);
                node.rho_fuel = rho_fuel(i,j,k);
                node.rho_co2 = rho_co2(i,j,k);
                node.rho_nox = rho_nox(i,j,k);
                node.rho_h2o = rho_h2o(i,j,k);
                node.temperature = temperature(i,j,k);
                node.pressure = calculate_pressure(&node);
                node.internal_energy = temperature_to_internal_energy(node.temperature);
                node.viscosity = temperature_to_viscosity(node.temperature);
                node.conductivity = temperature_to_conductivity(node.temperature);
                node.u = u(i,j,k);
                node.v = v(i,j,k);
                node.w = w(i,j,k);
                file << node.rho_o2 << " ";
                file << node.rho_n2 << " ";
                file << node.rho_fuel << " ";
                file << node.rho_co2 << " ";
                file << node.rho_nox << " ";
                file << node.rho_h2o << " ";
                file << node.pressure << " ";
                file << node.temperature << " ";
                file << node.viscosity << " ";
                file << node.internal_energy << " ";
                file << node.conductivity << " ";
                file << node.u << " ";
                file << node.v << " ";
                file << node.w << std::endl;
            }
        }
    }


  return 0;
}