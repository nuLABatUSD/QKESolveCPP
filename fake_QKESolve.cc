#include <iostream>
#include <cmath>
#include <chrono>
#include "mpi.h"
#include "globals.hh"
using namespace std::chrono;

#include "constants.hh"
#include "QKESolve.hh"
#include "QKE_methods.hh"

using std::cout;
using std::endl;

/**************************
/ Solving QKEs, just coherent terms
/ 
/ To compile:
/ mpic++ fake_QKESolve.cc QKE_methods.cc array_methods.cc thermodynamics.cc matrices.cc -std=c++11 -o qke
***************************/

//every processor meeds...
//integration object
//
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
    int sender, tag;
    MPI_Status* status;
    double myans;
    dummy_vars* E = d1->get_E();
    if(Globals::myid == 0){
        d2->zeros();
        
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
        
        
        
        for(int i=0; i<8*E->get_len(); i++){
            MPI_Recv(&myans, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status);
            sender = status->MPI_SOURCE;
            tag = status->MPI_TAG;
            d2->add_to_value(tag, myans);
            
        }
    }
    
    
    else{
        //what do i do about temperature? right now f doesn't do anything to it bc we assumed f to be 0 for temp
        //should i iterate through d1 instead of through eps? probably
        //i could change integrate to return a vector of all of p0 px py and pz for one p1 value
        //i could add a function so that given a certain index it can tell you if you have a neutrino or not and if it is p0, px, py, pz
        
        for(int i=0; i<E->get_len(); i+=Globals::numprocs){
            //neutrino
            //integration* new_integral = new integration(E, i);
            myans = 0;
            //myans = new_integral->whole_integral(d1, true, 0);
            MPI_Send(&myans, 1, MPI_DOUBLE, 0, 4*i, MPI_COMM_WORLD);
            //myans = new_integral->whole_integral(d1, true, 1);
            MPI_Send(&myans, 1, MPI_DOUBLE, 0, 4*i+1, MPI_COMM_WORLD);
            //myans = new_integral->whole_integral(d1, true, 2);
            MPI_Send(&myans, 1, MPI_DOUBLE, 0, 4*i+2, MPI_COMM_WORLD);
            //myans = new_integral->whole_integral(d1, true, 3);
            MPI_Send(&myans, 1, MPI_DOUBLE, 0, 4*i+3, MPI_COMM_WORLD);

            //antineutrino
            //myans = new_integral->whole_integral(d1, false, 0);
            MPI_Send(&myans, 1, MPI_DOUBLE, 0, 4*E->get_len()+4*i, MPI_COMM_WORLD);
            //myans = new_integral->whole_integral(d1, false, 1);
            MPI_Send(&myans, 1, MPI_DOUBLE, 0, 4*E->get_len()+4*i+1, MPI_COMM_WORLD);
            //myans = new_integral->whole_integral(d1, false, 2);
            MPI_Send(&myans, 1, MPI_DOUBLE, 0, 4*E->get_len()+4*i+2, MPI_COMM_WORLD);
            //myans = new_integral->whole_integral(d1, false, 3);
            MPI_Send(&myans, 1, MPI_DOUBLE, 0, 4*E->get_len()+4*i+3, MPI_COMM_WORLD);
            //delete new_integral;
        }
        
        
        
    }
}

int main(int argc, char *argv[])
{
    linspace_and_gl* et = new linspace_and_gl(0., 10., 201, 10);
    double eta_e = 0.01;
    double eta_mu = -0.01;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &(Globals::numprocs));
    MPI_Comm_rank(MPI_COMM_WORLD, &(Globals::myid));

    QKE* sim1 = new QKE(et, 0.8, 2.5e-15, eta_e, eta_mu);
    density* den1 = new density(et, eta_e, eta_mu);
    den1->set_T(0.25);
    sim1->set_ics(0, den1, 1.e12);
    sim1->run(1000, 10, 5.e15,"QKE1.csv", true);
    
    
    //run sim1 only in main but then f has to broadcast all necessary things
        
    
    MPI_Finalize();

    return 1;
        
}