#include "arrays.hh"
#include "QKE_methods.hh"

int main(){
    
    linspace_and_gl* eps = new linspace_and_gl(0,20,201,5);
    
    
    density* dens = new density(eps, 0.01, -0.01);
    dens->set_T(16.0);
    
    
    for(int i=1; i<50; i++){
        //if(i==26){
        nu_nu_annihilation* inte = new nu_nu_annihilation(eps, i, 32.0);
        //nu_e_collision* inte = new nu_e_collision(eps, i, 32.0);
        //inte->all_F_for_p1(dens,true, true);

        double* results = new double[4]();
        inte->whole_integral(dens, true, results, true);
        std::cout << results[0] << std::endl;

        delete inte;
        delete[] results;//}
        
    }
    
    
    
    
    delete eps;
    delete dens;
    
    return 0;
}