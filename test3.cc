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
    double eta_e = 0.01;
    double eta_mu = -0.01;
    
    
    QKESolveMPI* sim1 = new QKESolveMPI(myid, numprocs, et, 0.8, 2.5e-15, eta_e, eta_mu);
    density* den1 = new density(et, eta_e, eta_mu);
    density* den2 = new density(den1);
    den1->set_T(0.25);
    sim1->set_ics(0, den1, 1.e12);



    auto start = high_resolution_clock::now();
    sim1->f(1, den1, den2);

    //sim1->run(2, 1, 5.e15,"QKE1.csv", true);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    double time_elapsed = duration.count()/1000.;
    if(myid==0){
        den2->print_all();
        cout << "time elapsed: " << time_elapsed;
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

    //delete sim1;
    delete den2;
    delete den1;
    delete et;
    MPI_Finalize();
    //delete v;
    //delete test_int;
    //delete[] int_vals;
    
    return 0;
    
}