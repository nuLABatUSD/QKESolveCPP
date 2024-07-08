#include <iostream>
#include <cmath>
#include <chrono>
using namespace std::chrono;

#include "constants.hh"
#include "QKESolve.hh"
#include "QKE_methods.hh"

using std::cout;
using std::endl;

QKE::QKE(dummy_vars* epsilon, double sin_2theta, double delta_m_squared, double eta_e, double eta_mu) : ODESolve()
{
    sin_2theta = sin_2theta;
    cos_2theta = sqrt(1 - pow(sin_2theta, 2));
    delta_m_squared = delta_m_squared;

    y_values = new density(epsilon, eta_e, eta_mu);

    dummy_v_vac = new three_vector_for_QKE;
    dummy_v_vac->v_vacuum(delta_m_squared, cos_2theta, sin_2theta );

}

QKE::~QKE()
{ delete dummy_v_vac; }

void QKE::f(double t, density* d1, density* d2)
{
    d2->zeros();
    
    dummy_vars* E = d1->get_E();
    three_vector_for_QKE* dummy_v_dens = new three_vector_for_QKE;
    three_vector_for_QKE* dummy_v_therm = new three_vector_for_QKE;

    three_vector* V_nu = new three_vector;
    three_vector* V_nubar = new three_vector;
    three_vector* p = new three_vector;
    three_vector* vcrossp = new three_vector;
    
    dummy_v_dens->v_density(E, d1);
    dummy_v_therm->v_thermal(E, d1);
    double Tcm = d1->get_Tcm();
    
    double en = 0.;
    for (int i=1; i< E->get_len(); i++){
        en = E->get_value(i) * Tcm;
        
        V_nu->copy(dummy_v_dens);
        V_nu->add_to(1./en, dummy_v_vac);
        V_nu->add_to(en, dummy_v_therm);

        d1->p_vector(i, true, p);

        vcrossp->set_cross_product(V_nu, p);
        d2->set_value(4*i+1, vcrossp->get_value(0));
        d2->set_value(4*i+2, vcrossp->get_value(1));
        d2->set_value(4*i+3, vcrossp->get_value(2));

        V_nubar->copy(dummy_v_dens);
        V_nubar->add_to(-1./en, dummy_v_vac);
        V_nubar->add_to(-en, dummy_v_therm);

        d1->p_vector(i, false, p);
        vcrossp->set_cross_product(V_nubar, p);
        d2->set_value(4*(E->get_len())+4*i+1, vcrossp->get_value(0));
        d2->set_value(4*(E->get_len())+4*i+2, vcrossp->get_value(1));
        d2->set_value(4*(E->get_len())+4*i+3, vcrossp->get_value(2));

        

    }
        
    delete dummy_v_dens;
    delete dummy_v_therm;
    delete V_nu;
    delete V_nubar;
    delete vcrossp;
    delete p;
}

int main()
{    
    linspace_and_gl* et = new linspace_and_gl(0., 10., 201, 10);
    //linspace_for_trap* et = new linspace_for_trap(0., 10., 401);
    double eta_e = 0.01;
    double eta_mu = -0.01;
    /*
    QKE* sim1 = new QKE(et, 0.8, 2.5e-15, eta_e, eta_mu);
    density* den1 = new density(et, eta_e, eta_mu);

    den1->set_T(0.25);
    sim1->set_ics(0, den1, 1.e12);
    sim1->run(1000, 10, 5.e15,"QKE1-new.csv", true);

    QKE* sim2 = new QKE(et, 0.8, 2.5e-15, eta_e, eta_mu);
    density* den2 = new density(et, eta_e, eta_mu);

    den2->set_T(1.0);
    sim2->set_ics(0, den2, 1.e12);
    sim2->run(1000, 2000, 5.e18,"QKE2-new.csv", true);
    */
    QKE* sim3 = new QKE(et, 0.8, 2.5e-15, eta_e, eta_mu);
    density* den3 = new density(et, eta_e, eta_mu);

    den3->set_T(2.0);
    sim3->set_ics(0, den3, 1.e11);
    sim3->run(1000, 150000, 4.e19,"QKE3-new.csv", true);

    return 1;
        
}