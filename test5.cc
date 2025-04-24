#include "arrays.hh"
#include "QKE_methods.hh"

int main(){
    
    linspace_and_gl* eps = new linspace_and_gl(0,20,201,5);
    
    
    density* dens = new density(eps, 0.01, -0.01);
    dens->set_T(16.0);
    
    
    for(int i=0; i<eps->get_len(); i++){
        if(i==12){
        nu_e_collision* inte = new nu_e_collision(eps, i, 32.0);
        double* results = new double[4]();
        inte->whole_integral(dens, true, results, true);
        std::cout << results[3] << std::endl;
        delete inte;
        delete[] results;
        }
    }
    
    
    
    
    delete eps;
    delete dens;
    
    return 0;
}