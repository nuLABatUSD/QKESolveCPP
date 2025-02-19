#include "matrices.hh"
#include <iostream>
#include <fstream>
#include "arrays.hh"
#include <complex> 
#include <chrono>
#include "matrices.hh"
#include "constants.hh"
#include "QKE_methods.hh"
#include "alternative_integrals.hh"
#include <cmath>
#include <iomanip>

using std::cout;
using std::endl;
using std::complex;

int main(){
    

    
    linspace_and_gl_booles* bbb = new linspace_and_gl_booles(0,20,401,0);
    density* dens = new density(bbb, 0.01, -0.01);
    dens->set_T(1.0);
    for(int i=75; i<76; i++){
        nu_nu_collision* integral = new nu_nu_collision(bbb, i);
        double* results = new double[4]();

        integral->whole_integral(dens, true, results, true);
        std::cout << results[0] << std::endl;
        delete integral;
        delete[] results;
    }
    
    delete bbb;
    delete dens;
    
                                
    return 0;
}


