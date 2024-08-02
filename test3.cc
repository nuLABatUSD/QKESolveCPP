#include <iostream>
#include "arrays.hh"
#include "QKE_methods.hh"

using std::cout;
using std::endl;

int main(){
    linspace_and_gl* et = new linspace_and_gl(0., 10., 201, 5);
    double eta_e = 0.01;
    double eta_mu = -0.01;

    density* den1 = new density(et, eta_e, eta_mu);
    den1->set_T(0.25);
    density* den2 = new density(den1->num_bins(), et);
    double* int_vals = new double[4];

    for(int i=0; i<et->get_len(); i++){
        integration* test_int = new integration(et, i);
        test_int->whole_integral(den1, true, int_vals);
        for(int i=0; i<4; i++){
            cout << int_vals[i] << endl;
        }
        delete test_int;
    }
    
    

    
    /*
    test_int->whole_integral(den1, true, int_vals);
    cout << "neutrino results" << endl;
    for(int i=0; i<4; i++){
        cout << int_vals[i] << endl;
    }
    test_int->whole_integral(den1, false, int_vals);
    cout << "antineutrino results" << endl;
    for(int i=0; i<4; i++){
        cout << int_vals[i] << endl;
    }*/
/*
    three_vector* v = new three_vector();
    den1->p_vector(206, false, v);
    v->print_all();
   
    cout << den1->p0(205, false) << endl;

    
    double* dummy_int = new double[4];
    integration* test_int2 = new integration(et, 9);
    test_int2->whole_integral(den1, false, dummy_int);
    cout << "when processor 0 tries the integral" << endl;
    for(int j=0; j<4; j++){
        cout << dummy_int[j] << endl;
    }
    delete[] dummy_int;
    delete test_int2;*/


    delete den2;
    delete den1;
    delete et;
    //delete v;
    //delete test_int;
    delete[] int_vals;
    
    return 0;
    
}