#include <iostream>
#include "arrays.hh"
#include "QKE_methods.hh"

int main(){
    linspace_and_gl* et = new linspace_and_gl(0., 10., 201, 5);
    double eta_e = 0.01;
    double eta_mu = -0.01;
    
    density* den1 = new density(et, eta_e, eta_mu);
    density* den2 = new density(den1->num_bins(), et);
    integration* int0 = new integration(et, 1);
    double* int_vals = new double[4];
    int0->whole_integral(den1, true, int_vals);
    for(int i=0; i<4; i++){
        std::cout << int_vals[i] << std::endl;
    }
    
    delete den2;
    delete den1;
    delete et;
    delete int0;
    delete[] int_vals;
    
    return 0;
}