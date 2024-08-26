//#include "matrices.hh"
#include <iostream>
#include "arrays.hh"
#include <complex> 
#include "matrices.hh"
#include "constants.hh"
#include "QKE_methods.hh"

using std::cout;
using std::endl;
using std::complex;

int main(){
    int numbins = 201;
    linspace_and_gl* new_et = new linspace_and_gl(0,20,numbins,10);
    //std::cout << "numbins=" << numbins << ", p1=" << new_et->get_value(50) << std::endl;
    nu_e_collision* new_collision = new nu_e_collision(new_et, 0, .25);
    density* dens = new density(new_et, 0.01, -0.01);
    
    double* test_vals = new double[4]();
    new_collision->whole_integral(dens, true, test_vals);
    
    for(int i=0; i<4; i++){
        std::cout << "obtained: " << test_vals[i] << std::endl;}
    
    delete new_et;
    delete new_collision;
    delete dens;
    delete[] test_vals;
    
    return 0;
}
