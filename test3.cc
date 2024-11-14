#include <iostream>
#include "arrays.hh"
#include "QKE_methods.hh"
#include <chrono>
#include "alternative_integrals.hh"
#include <iomanip>
#include <fstream>

using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
    const std::string& term1_output = std::string(argv[1]);
    const std::string& term2_output = std::string(argv[2]);
    const std::string& p4_vals_output = std::string(argv[4]);
    const std::string& p3_vals_output = std::string(argv[3]);
    
    
    
    linspace_and_gl* eps = new linspace_and_gl(0,10,201,5);
    
    density* dens = new density(eps, 1, 8);
    dens->set_T(1);
    double* nu_nu_int = new double[4]();
    bool neutrino = true;
    
    int p1 = 75;
    int p2 = 205;
    double F0 = 0;
    three_vector* F = new three_vector();
    nu_nu_collision_one* og_collision = new nu_nu_collision_one(eps, p1);
    
    
    
    std::ofstream term1;
    term1.open(term1_output);

    std::ofstream term2;
    term2.open(term2_output);
    
    
    std::ofstream p4_vals;
    p4_vals.open(p4_vals_output);
    
    
    std::ofstream p3_vals;
    p3_vals.open(p3_vals_output);
    
    
    matrix* p_4 = new matrix();
    three_vector* dummy_p = new three_vector();
    
    double interpolated_p0 = 0;
    for(int p3=eps->get_len()-1; p3>=0; p3--){
       
        
        og_collision->Fvvsc_components_term_1(dens, neutrino, p2, p3, &F0, F);
        term1 << F->get_value(2) << ", ";
        
        og_collision->Fvvsc_components_term_2(dens, neutrino, p2, p3, &F0, F);
        term2 << F->get_value(2) << ", ";
        
        
        p3_vals << eps->get_value(p3) << ", ";
        
        
        double p4_energy = eps->get_value(p1)+ eps->get_value(p2) - eps->get_value(p3);
        p4_vals << p4_energy << ", ";
        
        /*
        
        p_4->convert_p4_to_interpolated_matrix(dens, neutrino, p4_energy);
        interpolated_p0 = real(p_4->get_A()->get_value(2));
        
        term1 << interpolated_p0 << ", ";
        dens->p0_p(p3, neutrino, dummy_p);
        term2 << dummy_p->get_value(2) << ", ";
        
        */
        
    }
    
    
    term1.close();
    term2.close();
    p4_vals.close();
    p3_vals.close();
        
    delete F;
    delete og_collision;
    delete p_4;
    delete dummy_p;
    
    delete eps;
    delete dens;
    delete[] nu_nu_int;
    return 0;
    
}