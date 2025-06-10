#include <cmath>
#include <iostream>
#include <chrono>
#include <fstream>
#include <ostream>

#include "constants.hh"
#include "arrays.hh"
#include "QKE_methods.hh"

using namespace std;

int main(int argc, char* argv[]){
    const std::string& output_file = std::string(argv[1]);
    
    std::ofstream output;
    output.open(output_file);
    
    linspace_and_gl* eps = new linspace_and_gl(0, 20, 201, 5);
    density* dens = new density(eps, 0.01, -0.01);
    dens->set_T(1.0);
    double* net_results = new double[4]();
    double* FRS_results = new double[4]();
    
    bool neutrino = true;
    double ratio;
    
    for (int p1=0; p1<eps->get_len(); p1++){
        auto start = std::chrono::high_resolution_clock::now();
        nu_nu_annihilation* integral = new nu_nu_annihilation(eps, p1, 32);
        integral->whole_integral(dens, neutrino, net_results, true);
        integral->whole_integral(dens, neutrino, FRS_results, false);
        
        if (FRS_results[0]==0){
            ratio = 0;
        }
        else{
            ratio = net_results[0] / FRS_results[0];
        }
        auto stop = std::chrono::high_resolution_clock::now();
        
        auto duration = duration_cast<std::chrono::milliseconds>(stop - start);
        double time_elapsed = duration.count()/1000.;
            
        std::cout << "p1 energy is " << eps->get_value(p1) << ": " << ratio << ", found in " << time_elapsed << " seconds." << std::endl;
        output << ratio << ", ";
        delete integral;
        
    }
    output.close();
    
    delete eps;
    delete dens;
    delete[] FRS_results;
    delete[] net_results;
    
    
    return 0;
}