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
    
    linspace_and_gl* eps = new linspace_and_gl(0,20,201,5);
    density* dens = new density(eps, 0.01, -0.01);
    dens->set_T(1.0);
    
    nu_e_collision* col = new nu_e_collision(eps, 100, 32.0);
    std::cout << "net_true = [";
    col->all_F_for_p1(dens, true, true);
    std::cout << "]" << std::endl;
    /*
    std::cout << "net_false = [";
    col->all_F_for_p1(dens, true, false);
    std::cout << "]" << std::endl;
    */
    
    
    delete eps;
    delete dens;
    delete col;
                     
    return 0;
}


