#include "arrays.hh"
#include "QKE_methods.hh"
#include <chrono>

int main(){
    
    linspace_and_gl* eps = new linspace_and_gl(0,20,201,5);
    linspace_and_gl* eps_finer = new linspace_and_gl(0,40,4001,0);
    density* dens = new density(eps, 0.01, -0.01);
    density* dens_finer = new density(eps_finer, 0.01, -0.01);
    three_vector* p = new three_vector();
    
    
    double err;
    for(int i=0; i<eps_finer->get_len(); i++){
        err = pow(dens->interpolated_matrix(true, eps_finer->get_value(i), p)-dens_finer->p0(i, true), 2);
        std::cout << err << ", ";
    }
    
    
    delete eps;
    delete dens;
    delete eps_finer;
    delete dens_finer;
    delete p;
    return 0;
}