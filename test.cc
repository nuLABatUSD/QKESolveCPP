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

using std::cout;
using std::endl;
using std::complex;

int main(){
    int numbins = 201;
    linspace_and_gl* new_et = new linspace_and_gl(0,10,numbins,5);
    density* dens = new density(new_et, 0.01, -0.01);
    dens->set_T(1);
    bool neutrino = true;
    
    
    // print R 
    
    // integral_one is regular (-)
    // integral_one_1 is plus (+)
    
    int p1 = 75;
    //int p2 = 50;
    
    for(int p2=0; p2<new_et->get_len(); p2++){
        nu_nu_collision_one* regular = new nu_nu_collision_one(new_et, p1);
        regular->Fvvsc_for_p1(dens, neutrino);
        regular->Fvvbarsc_for_p1(dens, neutrino);
        double* regular_results = new double[4]();
        double int_one = regular->interior_integral(p2, 0);
        std::cout << "done w integral one" << std::endl;
    
        nu_nu_collision_one_1* plus = new nu_nu_collision_one_1(new_et, p1);
        plus->Fvvsc_for_p1(dens, neutrino);
        plus->Fvvbarsc_for_p1(dens, neutrino);
        double* plus_results = new double[4]();
        double int_two = plus->interior_integral(p2,0);
        std::cout << "done w integral two" << std::endl;
        
        std::cout << int_one / int_two << std::endl;
    
        //std::cout << regular_results[0] / plus_results[0] << std::endl;
        delete[] regular_results;
        delete[] plus_results;
        delete regular;
        delete plus;
        
    }
    
    
    
    
    
    delete dens;
    delete new_et;
    
    return 0;
}


