//#include "matrices.hh"
#include <iostream>
#include <fstream>
#include "arrays.hh"
#include <complex> 
#include "matrices.hh"
#include "constants.hh"
#include <iomanip>
#include "gl_vals.hh"
#include <string>

using std::cout;
using std::endl;
using std::complex;

int main(){
    
    
    linspace_for_trap* et = new linspace_for_trap(0.,20, 201);
    double eta_e = 0.2;
    double eta_mu = -0.02;
    linspace_and_gl* new_et = new linspace_and_gl(0,20,201,5);
    
    dep_vars* dens_vals = new dep_vars(8*new_et->get_len()+2);
    
    std::ifstream densfile;
    densfile.open("output2.csv");
    if (!densfile.is_open()) {
        std::cout << "Error opening density input file" << std::endl;
    }
    
    std::string mystring;
    //for(int i=0; i<dens_vals->length(); i++){
    int i=0;
    while(densfile){
        std::getline(densfile, mystring, ',');
        dens_vals->set_value(i, std::stod(mystring));
        i++;
    }
    densfile.close();
    for(int i=0; i<dens_vals->length(); i++){
        cout << dens_vals->get_value(i) << endl;
    }
    
    delete et;
    delete new_et;
    //delete dens_vals;
    
    return 0;
}
