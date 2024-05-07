#include <iostream>
#include <cmath>
#include "quad_vals.hh"
#include "constants.hh"
#include <fstream>
#include "energy_density_and_pressure.hh"

using std::cout;
using std::endl;
using std::ofstream;

int main(){
    
    ofstream file("energy_density_and_pressure_for_T.csv");
    
    file << "T, density, pressure" << endl;
    
    for (int i=1; i<101; i++){
        double dens = 0.;
        double press =0.;
        energy_and_pressure(_electron_mass_, i, &dens, &press);
        
        file << i << ", " << dens << ", " << press << endl;
    }
    file.close();
    
}
