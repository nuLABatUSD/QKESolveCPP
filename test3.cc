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
    const std::string& initialdens = std::string(argv[5]);
    
    
    
    linspace_and_gl* eps = new linspace_and_gl(0,20,201,0);
    density* dens = new density(eps, 0.01, -0.01);
    /*
    for(int i=0; i<dens->length(); i++){
        dens->set_value(i, dens->get_value(i));
        
    }
    */
    /*
    double* dens_vals = new double[8*eps->get_len()+2]();
    
    std::ifstream densfile;
    densfile.open(initialdens);
    if (!densfile.is_open()) {
        std::cout << "Error opening density input file" << std::endl;
    }
    std::string mystring;
    int i=0;
    while(densfile){
        std::getline(densfile, mystring, ',');
        dens_vals[i] = std::stod(mystring);
        i++;
    }
    densfile.close();
    density* dens = new density(eps->get_len(), eps, dens_vals);
    for(int i=0; i<dens->length(); i++){
        dens->set_value(i, dens->get_value(i)/10);
        
    }
    delete[] dens_vals;
    */
    dens->set_T(1);
    double* nu_nu_int = new double[4]();
    bool neutrino = true;
    
    int p1 = 50;
    //linspace_and_gel* eeps = new linspace_and_gel(eps, eps->get_value(p2)+eps->get_value(p1), 10);
    double F0 = 0;
    three_vector* F = new three_vector();
    nu_nu_collision* new_collision = new nu_nu_collision(eps, p1);
    nu_nu_collision* og_collision = new nu_nu_collision(eps, p1);
    nu_nu_collision_one* new_collision_1 = new nu_nu_collision_one(eps, p1);
    std::cout << "p1 energy=" << eps->get_value(p1) << std::endl;
    
    
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
    
    double* int_vals = new double[4]();
    
    double interpolated_p0 = 0;
    
    
    new_collision->Fvvsc_for_p1(dens, neutrino);
    new_collision->Fvvbarsc_for_p1(dens, neutrino);
    new_collision_1->Fvvsc_for_p1(dens, neutrino);
    new_collision_1->Fvvbarsc_for_p1(dens, neutrino);
    
    
    /*
    for(int p2=0; p2<eps->get_len(); p2++){
        new_collision->interior_integral(p2, 0);
    }*/
    /*
    std::cout << "------------" << std::endl;
    for(int p2=0; p2<eps->get_len(); p2++){
        std::cout << eps->get_value(p2) << ", ";
    }*/
    int p2=10;
    //for(int p3=0; p3<eps->get_len(); p3++){
    for(int i=0; i<eps->get_len(); i++){
        /*
        double p4_energy = eps->get_value(p1) + eps->get_value(p2) - eps->get_value(p3);
        
        
        if (p4_energy >=0){
            
            new_collision_1->Fvvbarsc_components(dens, neutrino, p2, p3, &F0, F);
            term1 << F0 << ", ";

            new_collision_1->Fvvbarsc_components_term_2(dens, neutrino, p2, p3, &F0, F);
            term2 << F0 << ", ";
            p3_vals << eps->get_value(p3) << ", ";
            p4_vals << p4_energy << ", ";
        }
        
        */

            term1 << new_collision->interior_integral(i, 0) << ", ";
            term2 << new_collision_1->interior_integral(i, 0) << ", ";
            
            p4_vals << eps->get_value(i) << ", ";
            p3_vals << eps->get_value(i) << ", ";
            //p3_vals << eps->get_value(p1) + eps->get_value(p2) - eeps->get_value(p3) << ", ";

            /*
            p4_energy = eps->get_value(p1)+ eps->get_value(p2) - eps->get_value(p3);
            p4_vals << p4_energy << ", ";
            p3_vals << eps->get_value(p3) << ", ";
            */
       // }
        /*
        
        p3_vals << eps->get_value(p3) << ", ";
        p4_vals << p4_energy << ", ";
        
        
        interpolated_p0 = dens->interpolate_p0(neutrino, p4_energy);
        term1 << interpolated_p0 << ", ";
        dens->p0_p(p3, neutrino, dummy_p);
        term2 << dens->p0(p3, neutrino) << ", ";
        
        */
        
    }

    /*
    og_collision->Fvvsc_for_p1(dens, neutrino);
    og_collision->Fvvbarsc_for_p1(dens, neutrino);
    
    new_collision->Fvvsc_for_p1(dens, neutrino);
    new_collision->Fvvbarsc_for_p1(dens, neutrino);
    
    for(int p2=0; p2<eps->get_len(); p2++){
        p3_vals << eps->get_value(p2) << ", ";
        p4_vals << eps->get_value(p2) << ", ";
        
        term1 << og_collision->interior_integral(p2, 0) << ", ";
        
        term2 << new_collision->interior_integral(p2, 0) << ", ";
        
        
    }*/
    delete[] int_vals;
    
    term1.close();
    term2.close();
    p4_vals.close();
    p3_vals.close();
        
    delete F;
    delete og_collision;
    delete new_collision;
    delete new_collision_1;
    //delete p_4;
    //delete dummy_p;
    
    delete eps;
    //delete eeps;
    delete dens;
    delete[] nu_nu_int;
    return 0;
    
}