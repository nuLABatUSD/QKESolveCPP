#include <iostream>
#include "QKESolveMPI.hh"
#include "arrays.hh"
#include "QKE_methods.hh"

#include "mpi.h"
using std::cout;
using std::endl;
using namespace std;


/*
TO RUN:
mpic++ test2.cc QKESolveMPI.cc array_methods.cc QKESolve.cc QKE_methods.cc thermodynamics.cc matrices.cc -std=c++11 -o wed
mpiexec -n 4 wed

*/

int main(int argc, char *argv[])
{
    int myid, numprocs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    linspace_and_gl* et = new linspace_and_gl(0., 10., 201, 5);
    //linspace_for_trap* et = new linspace_for_trap(0., 10., 401);
    double eta_e = 0.01;
    double eta_mu = -0.01;
    
    QKESolveMPI* sim1 = new QKESolveMPI(myid, numprocs, et, 0.8, 2.5e-15, eta_e, eta_mu);
    density* den1 = new density(et, eta_e, eta_mu);
    density* den2 = new density(den1->num_bins(), et);
    den1->set_T(0.25);
    sim1->set_ics(0, den1, 1.e12);
    sim1->f(1, den1, den2);
    
    MPI_Finalize();
    
    delete et;
    delete sim1;
    delete den1;
    delete den2;
    return 0;
    
}

