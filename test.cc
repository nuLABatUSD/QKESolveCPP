#include "matrices.hh"
#include <iostream>
#include <fstream>
#include "arrays.hh"
#include <complex> 
#include <chrono>
#include "matrices.hh"
#include "constants.hh"
#include "QKE_methods.hh"
#include "thermodynamics.hh"
#include "alternative_integrals.hh"
#include <cmath>
#include <iomanip>

using std::cout;
using std::endl;
using std::complex;


int main(int argc, char* argv[]){
    const std::string& input_file = std::string(argv[1]);
    
    linspace_and_gl* epsilon = new linspace_and_gl(0,20,201,5);
    double* dens_vals = new double[epsilon->get_len()*8+2]();
    
    std::ifstream densfile;
    densfile.open(input_file);


    if (!densfile.is_open()) {
        std::cout << "Error opening density input file" << std::endl;
    }
    
    
    int j = 0;
    std::string line;
    while(std::getline(densfile, line)){
        std::string densval;
        std::string delimiter = ", ";

        size_t pos = 0;
        int i=0;
        while((pos = line.find(delimiter)) != std::string::npos){
            densval = line.substr(0, pos);
            //this takes care of first two elements being initial place and initial step
            if(i>1){
                dens_vals[i-2] = std::stod(densval);
            }
            line.erase(0, pos + delimiter.length());
            i++;
        }
        dens_vals[i-2] = std::stod(line);
        density* dens = new density(epsilon->get_len(), epsilon, dens_vals);
        
        double* endens = new double[4]();
        dens->energy_density(endens);
        double nd=0;
        for(int i=0; i<4; i++){
            nd += endens[i];
        }
        delete[] endens;
        
        double rho;
        double p;
        energy_and_pressure(_electron_mass_, dens->get_Tcm(), &rho, &p);
        
        std::cout << "On line " << j << ", neutrino energy density is " << nd << " and electron energy density is " << rho << ", sum is " << rho+nd << std::endl;
        
        delete dens;

        j++;

    }
    densfile.close();
    
    
    delete[] dens_vals;
    delete epsilon;
    
    
    return 0; 
}
