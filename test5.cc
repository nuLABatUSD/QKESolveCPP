#include "arrays.hh"
#include "QKE_methods.hh"

int main(){
    
    linspace_and_gl* eps = new linspace_and_gl(0,20,201,5);
    density* dens = new density(eps, 1, 8);
    dens->set_T(1.0);
    
    nu_nu_collision* inte = new nu_nu_collision(eps, 201);
    double* results = new double[4]();
    inte->whole_integral(dens, true, results, true);
    std::cout << results[0] << std::endl;
    
    delete eps;
    delete dens;
    delete inte;
    delete[] results;
    return 0;
}