#include <cmath>
#include "quad_vals.hh"
#include "constants.hh"

double energy_f(double u, double m, double T){
    return 2 / (pow(_PI_,2)) * pow(u,2) * pow(T,3) * sqrt(pow(u,2)*pow(T,2) + pow(m,2)) / (exp(sqrt(pow(u,2)*pow(T,2) + pow(m,2))/T) + 1) * exp(u);
}

double pressure_f(double u, double m, double T){
    return 2 / (3 * pow(_PI_,2)) * pow(u,4) * pow(T,5) / sqrt(pow(u,2)*pow(T,2) + pow(m,2)) / (exp(sqrt(pow(u,2)*pow(T,2) + pow(m,2))/T) + 1) * exp(u);
}
                                                               
void energy_and_pressure(double m, double T, double* rho, double* P){
    *rho = 0;
    *P = 0;
    for (int i=0; i<nvals; i++){
        *rho += energy_f(xvals[i],m,T) * wvals[i];
        *P += pressure_f(xvals[i],m,T) * wvals[i];
    }
}