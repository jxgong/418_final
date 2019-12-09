#include "physics.h"
#include "common.h"
#include <math.h>

float calculate_pressure(Node *node){
    float nV = 0.f;
    nV += node->rho_o2 / o2_molar_mass;
    nV += node->rho_n2 / n2_molar_mass;
    nV += node->rho_fuel / fuel_molar_mass;
    nV += node->rho_co2 / co2_molar_mass;
    nV += node->rho_nox / nox_molar_mass;
    nV += node->rho_h2o / h2o_molar_mass;
    return nV * 8.31 * node->temperature;
}

float temperature_to_viscosity(float temp){
    return 0.000001458f * sqrtf(temp*temp*temp) / (temp + 110.4f);
}

float temperature_to_conductivity(float temp){
    return 35.f;
}

float internal_energy_to_temperature(float temperature){
    // TODO: enter empirical mapping of internal energy to temperature
    return boundary_idling_temp;
}

float temperature_to_internal_energy(float temperature){
    // TODO: enter empirical mapping of internal energy to temperature
    return 268.f;
}

float o2_enthalpy(float temperature){
    float T = temperature;
    if (temperature < 700.f){
        return 1000.f*(31.332f*T-20.235f*T*T+57.866f*T*T*T-36.506f*T*T*T*T-0.007f/T-8.903f-0.f);
    }
    else if (temperature < 2000.f){
        return 1000.f*(30.032f*T+8.773f*T*T-3.998f*T*T*T+0.788f*T*T*T*T-0.742f/T-11.325f-0.f);
    }
    else{
        return 1000.f*(20.911f*T+10.721f*T*T-2.020f*T*T*T+0.146f*T*T*T*T+9.246f/T+5.338f-0.f);
    }
}