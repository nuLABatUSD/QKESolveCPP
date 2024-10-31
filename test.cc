//#include "matrices.hh"
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

int main(int argc, char* argv[]){
    
    const std::string& og_output = std::string(argv[1]);
    const std::string& output_one = std::string(argv[2]);
    const std::string& output_two = std::string(argv[3]);
    
    
    int numbins = 201;
    linspace_and_gl* new_et = new linspace_and_gl(0,10,201,5);
    density* dens = new density(new_et, 0.01, -0.01);
    dens->set_T(1);
    bool neutrino = true;
    
    dep_vars* og_int_vals = new dep_vars(new_et->get_len());
    dep_vars* int_vals_one = new dep_vars(new_et->get_len());
    dep_vars* int_vals_two = new dep_vars(new_et->get_len());
    
    double* nu_nu_int = new double[4]();
    
    
    std::ofstream ogfile;
    ogfile.open(og_output);

    std::ofstream fileone;
    fileone.open(output_one);

    std::ofstream filetwo;
    filetwo.open(output_two);
    
    for(int p1=0; p1<new_et->get_len(); p1++){
        auto start = std::chrono::high_resolution_clock::now();
        nu_nu_collision* og_collision = new nu_nu_collision(new_et, p1);
        og_collision->whole_integral(dens, neutrino, nu_nu_int);
        og_int_vals->set_value(p1, pow(new_et->get_value(p1),2) * nu_nu_int[0]);
        ogfile << og_int_vals->get_value(p1) << ", ";
        
        nu_nu_collision_one* collision_one = new nu_nu_collision_one(new_et, p1);
        collision_one->whole_integral(dens, neutrino, nu_nu_int);
        int_vals_one->set_value(p1, pow(new_et->get_value(p1),2) * nu_nu_int[0]);
        fileone << int_vals_one->get_value(p1) << ", ";
        
        nu_nu_collision_two* collision_two = new nu_nu_collision_two(new_et, p1);
        collision_two->whole_integral(dens, neutrino, nu_nu_int);
        int_vals_two->set_value(p1, pow(new_et->get_value(p1),2) * nu_nu_int[0]);
        filetwo << int_vals_two->get_value(p1) << ", ";
        
        delete og_collision;
        delete collision_one;
        delete collision_two;
        
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = duration_cast<std::chrono::milliseconds>(stop - start);
        double time_elapsed = duration.count()/1000.;
        
        std::cout << "Finished p1=" << new_et->get_value(p1) << " in " << time_elapsed << " seconds" << std::endl;
    }
    
    ogfile.close();
    fileone.close();
    filetwo.close();
    
    double og_result = new_et->integrate(og_int_vals);
    double result_one = new_et->integrate(int_vals_one);
    double result_two = new_et->integrate(int_vals_two);
    
    std::cout << "dn/dt for original method is " << og_result << std::endl;
    std::cout << "dn/dt for first new method is " << result_one << std::endl;
    std::cout << "dn/dt for second new method is " << result_two << std::endl;
    
    
    delete[] nu_nu_int;
    delete og_int_vals;
    delete int_vals_one;
    delete int_vals_two;
    delete dens;
    delete new_et;
    
    return 0;
}


