#include "arrays.hh"
#include "QKE_methods.hh"

int main(){
    
    linspace_and_gl* eps = new linspace_and_gl(0,20,201,5);
    density* dens = new density(eps, 1, 8);
    dens->set_T(16.0);
    
    for(int i=0; i<eps->get_len(); i++){
        nu_nu_collision* inte = new nu_nu_collision(eps, i);
        double* results = new double[4]();
        inte->whole_integral(dens, false, results, true);
        std::cout << results[3] << std::endl;
        delete inte;
        delete[] results;
    }
    
    
    
    
    delete eps;
    delete dens;
    
    return 0;
}