#include <iostream>
#include "arrays.hh"
#include "QKE_methods.hh"
#include <chrono>
#include "QKESolveMPI.hh"
#include "mpi.h"

using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
    int myid, numprocs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    linspace_and_gl* et = new linspace_and_gl(0., 10., 201, 5);
    double eta_e = -0.1;
    double eta_mu = 0.1;
    
    QKESolveMPI* sim1 = new QKESolveMPI(myid, numprocs, et, 0.8, 2.5e-15, eta_e, eta_mu);
    density* den1 = new density(et, eta_e, eta_mu);
    density* den2 = new density(den1);
    den1->set_T(0.25);
    sim1->set_ics(0, den1, 1.e12);
    
    sim1->f(1, den1, den2);
    
    if(myid==0){
        den2->print_all();
    }
    
    
    /*
    double* int_vals = new double[4]();
    integration* test_int = new integration(et, 1);
    test_int->whole_integral(den1, true, int_vals);
    if(myid==0){
        cout << "neutrino results" << endl;
        for(int i=0; i<4; i++){
            cout << int_vals[i] << endl;
        }
    }
    test_int->whole_integral(den1, false, int_vals);
    if(myid==0){
        cout << "antineutrino results" << endl;
        for(int i=0; i<4; i++){
            cout << int_vals[i] << endl;
        }
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

    delete sim1;
    delete den2;
    delete den1;
    delete et;
    MPI_Finalize();
    //delete v;
    //delete test_int;
    //delete[] int_vals;
    
    return 0;
    
}