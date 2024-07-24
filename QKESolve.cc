#include <iostream>
#include <cmath>
#include <chrono>
using namespace std::chrono;

#include "constants.hh"
#include "QKESolve.hh"
#include "QKE_methods.hh"

using std::cout;
using std::endl;

/**************************
/ Solving QKEs, just coherent terms -> now with collisions!
/ 
/ To compile:
/ g++ QKESolve.cc QKE_methods.cc array_methods.cc thermodynamics.cc matrices.cc -std=c++11 -o qke
***************************/
QKE::QKE(linspace_and_gl* e, double sin_2theta, double delta_m_squared, double eta_e, double eta_mu) : ODESolve()
{
    epsilon = new linspace_and_gl(e);
    sin_2theta = sin_2theta;
    cos_2theta = sqrt(1 - pow(sin_2theta, 2));
    delta_m_squared = delta_m_squared;

    //y_values = new density(epsilon, eta_e, eta_mu);

    dummy_v_vac = new three_vector_for_QKE;
    dummy_v_vac->v_vacuum(delta_m_squared, cos_2theta, sin_2theta);
    /*
    int_objects = new integration*[epsilon->get_len()];
    
    
    for(int i=0; i<epsilon->get_len(); i++){
        int_objects[i] = new integration(epsilon, i);
    }
    */

}

QKE::~QKE()
{
    delete dummy_v_vac;
    /*
    for(int i=0; i<epsilon->get_len(); i++){
        delete int_objects[i];
    }
    delete[] int_objects;*/
    cout << "the issue is not deleting dummy_v_vac" << endl;
    delete epsilon;

    


}


void QKE::f(double t, density* d1, density* d2)
{
    d2->zeros();
    
    //is E ever going to be different from epsilon
    three_vector_for_QKE* dummy_v_dens = new three_vector_for_QKE;
    three_vector_for_QKE* dummy_v_therm = new three_vector_for_QKE;

    three_vector* V_nu = new three_vector;
    three_vector* V_nubar = new three_vector;
    three_vector* p = new three_vector;
    three_vector* vcrossp = new three_vector;
    
    cout << "i got to here in f" << endl;
    dummy_v_dens->v_density(epsilon, d1);
    dummy_v_therm->v_thermal(epsilon, d1);
    double Tcm = d1->get_Tcm();
    
    double en = 0.;
    double* integral_vals = new double[4];
    for (int i=1; i< epsilon->get_len(); i++){
        en = epsilon->get_value(i) * Tcm;
        
        V_nu->copy(dummy_v_dens);
        V_nu->add_to(1./en, dummy_v_vac);
        V_nu->add_to(en, dummy_v_therm);

        d1->p_vector(i, true, p);

        vcrossp->set_cross_product(V_nu, p);
        d2->set_value(4*i+1, vcrossp->get_value(0));
        d2->set_value(4*i+2, vcrossp->get_value(1));
        d2->set_value(4*i+3, vcrossp->get_value(2));
        
        /*
        int_objects[i]->whole_integral(d1, true, integral_vals);
        for(int j=0; j<4; j++){
            d2->add_to_value(4*i+j, integral_vals[j]);
        }*/
        /*
        
        d2->set_value(4*i, int_objects[i]->whole_integral(d1, true, 0));
        d2->add_to_value(4*i+1, int_objects[i]->whole_integral(d1, true, 1));
        d2->add_to_value(4*i+2, int_objects[i]->whole_integral(d1, true, 2));
        d2->add_to_value(4*i+3, int_objects[i]->whole_integral(d1, true, 3));*/
        

        V_nubar->copy(dummy_v_dens);
        V_nubar->add_to(-1./en, dummy_v_vac);
        V_nubar->add_to(-en, dummy_v_therm);

        d1->p_vector(i, false, p);
        vcrossp->set_cross_product(V_nubar, p);
        d2->set_value(4*(epsilon->get_len())+4*i+1, vcrossp->get_value(0));
        d2->set_value(4*(epsilon->get_len())+4*i+2, vcrossp->get_value(1));
        d2->set_value(4*(epsilon->get_len())+4*i+3, vcrossp->get_value(2));
        /*
        int_objects[i]->whole_integral(d1, false, integral_vals);
        for(int j=0; j<4; j++){
            d2->add_to_value(4*epsilon->get_len()+4*i+j, integral_vals[j]);
        }*/

        /*
        d2->set_value(4*(E->get_len())+4*i, int_objects[i]->whole_integral(d1, false, 0));
        d2->add_to_value(4*(E->get_len())+4*i+1, int_objects[i]->whole_integral(d1, false, 1));
        d2->add_to_value(4*(E->get_len())+4*i+2, int_objects[i]->whole_integral(d1, false, 2));
        d2->add_to_value(4*(E->get_len())+4*i+3, int_objects[i]->whole_integral(d1, false, 3));*/

    }
    
    
    
    
    
        
    delete dummy_v_dens;
    delete dummy_v_therm;
    delete V_nu;
    delete V_nubar;
    delete vcrossp;
    delete p;
    delete[] integral_vals;
    
   
   
}

int main()
{    
    linspace_and_gl* et = new linspace_and_gl(0., 10., 201, 5);
    //linspace_for_trap* et = new linspace_for_trap(0., 10., 401);
    double eta_e = 0.01;
    double eta_mu = -0.01;
    
    // 4 seconds
    cout << "i got to before declaration of sim1" << endl; //yes
    QKE* sim1 = new QKE(et, 0.8, 2.5e-15, eta_e, eta_mu);
    cout << "i got past defining sim1" <<endl; //yes
    density* den1 = new density(et, eta_e, eta_mu);
    density* den2 = new density(den1->num_bins(), et);
    
    den1->set_T(0.25);
    sim1->set_ics(0, den1, 1.e12);
    //auto start = high_resolution_clock::now();
    cout << "i got to before f" << endl;
    sim1->f(1, den1, den2);
    //sim1->run(10, 1, 5.e15,"QKE1.csv", true);
    
    /*
    integration* test_int = new integration(et, 200);
    double* test = new double[4];
    test_int->whole_integral(den1, true, test);
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << endl << "Time elapsed: "
         << duration.count()/1000. << " seconds" << endl;
    */
    delete et;
    delete sim1;
    delete den1;
    delete den2; 
    //delete test_int;
    //delete[] test;

    /*
    // 21:44 (mm:ss)
    QKE* sim2 = new QKE(et, 0.8, 2.5e-15, eta_e, eta_mu);
    density* den2 = new density(et, eta_e, eta_mu);
    den2->set_T(1.0);
    sim2->set_ics(0, den2, 1.e12);
    sim2->run(1000, 2000, 5.e18,"QKE2.csv", true);

    
    // 22:01 (hh:mm)
    QKE* sim3 = new QKE(et, 0.8, 2.5e-15, eta_e, eta_mu);
    density* den3 = new density(et, eta_e, eta_mu);
    den3->set_T(2.0);
    sim3->set_ics(0, den3, 1.e11);
    sim3->run(1000, 150000, 4.e19,"QKE3.csv", true);
    */
    return 1;
        
}