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
    linspace_and_gl* new_et = new linspace_and_gl(0,10,201,5);
    
    for(int p1=0; p1<204; p1++){
        nu_e_collision* new_collision = new nu_e_collision(new_et, p1, 1);
        /*
        density* dens = new density(new_et, 0.01, -0.01);
        dens->set_T(1);

        double* test_vals = new double[4]();
        new_collision->whole_integral(dens, true, test_vals);

        std::cout << "p1=" << p1 << ", p0 int=" << test_vals[0] << std::endl;
        
        delete dens;*/
        delete new_collision;
        //delete[] test_vals;
            
    }
    
    delete new_et;
    
    return 0;
}
