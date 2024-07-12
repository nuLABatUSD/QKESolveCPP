#include <iostream>
#include "constants.hh"
#include "arrays.hh"
#include "QKE_methods.hh"

using std::cout;
using std::endl;

int main(){
    
    double eta_e = 0.2;
    double eta_mu = -0.02;
    linspace_and_gl* new_et = new linspace_and_gl(0,20,201,2);
    
    
    density* new_den = new density(new_et, eta_e, eta_mu);
    integration* iiii = new integration(new_et, 2);
    cout << iiii->whole_integral(new_den, true, 0);
    delete iiii;
    delete new_et;
    delete new_den;
    
    return 0;
}
