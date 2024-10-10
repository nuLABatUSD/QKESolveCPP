//#include "matrices.hh"
#include <iostream>
#include "arrays.hh"
#include <complex> 
#include "matrices.hh"
#include "constants.hh"
#include "QKE_methods.hh"
#include "alternative_integrals.hh"

using std::cout;
using std::endl;
using std::complex;

int main(){
    int numbins = 201;
    linspace_and_gl* new_et = new linspace_and_gl(0,10,201,5);
    
    //for(int p1=0; p1<204; p1++){
        nu_nu_collision_two* new_collision = new nu_nu_collision_two(new_et, 40);
        
        density* dens = new density(new_et, 0.01, -0.01);
        dens->set_T(1);

        double* test_vals = new double[4]();
        new_collision->whole_integral(dens, true, test_vals);

        
        for(int j=0; j<4; j++){
            std::cout << test_vals[j] << std::endl;
        }
        //std::cout << "p1=" << p1 << ", p0 int=" << test_vals[0] << std::endl;
        
        delete dens;
        delete new_collision;
        delete[] test_vals;
            
    //}
    
    delete new_et;
    
    return 0;
}
