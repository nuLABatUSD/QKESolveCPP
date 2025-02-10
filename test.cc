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
    
    linspace_and_gl_booles* bbb = new linspace_and_gl_booles(0,5,101,0);
    dep_vars* squares = new dep_vars(101);
    for (int i=0; i< squares->length(); i++){
        squares->set_value(i, pow(bbb->get_value(i),5));
    }
    
    std::cout << std::setprecision(16) << bbb->integrate(squares);
    
    return 0;
}


