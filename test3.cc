#include <iostream>
#include "arrays.hh"
#include "QKE_methods.hh"

int main(){
    linspace_and_gl* et = new linspace_and_gl(0., 10., 201, 5);
    double eta_e = 0.01;
    double eta_mu = -0.01;
    
    density* den1 = new density(et, eta_e, eta_mu);
    density* den2 = new density(den1->num_bins(), et);
    
    delete den2;
    delete den1;
    delete et;
    
    return 0;
}