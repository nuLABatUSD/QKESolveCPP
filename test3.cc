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
    double eta_e = -0.01;
    double eta_mu = 0.01;
    
    QKESolveMPI* sim1 = new QKESolveMPI(myid, numprocs, et, 0.8, 2.5e-15, 0, 1.e12, "initialthermaldistribution.csv");
    density* den1 = new density(et, eta_e, eta_mu);
    density* den2 = new density(den1);
    den1->set_T(0.25);
    sim1->set_ics(0, den1, 1.e12);
    
    double average_time_elapsed = 0;
    for(int i=0; i<5; i++){
        auto start = high_resolution_clock::now();
        sim1->f(1, den1, den2);
        auto stop = high_resolution_clock::now();
        
        auto duration = duration_cast<milliseconds>(stop - start);
        double time_elapsed = duration.count()/1000.;
        double max_time_elapsed = 0;
        MPI_Reduce(&time_elapsed, &max_time_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if(myid==0){std::cout << "Time for trial " << i << ": " << max_time_elapsed << std::endl;}
        average_time_elapsed += max_time_elapsed;
    }
    average_time_elapsed /= 5;
    if(myid==0){
        std::cout << "Average Time Elapsed in 5 Trials: " << average_time_elapsed << std::endl;
    }
    delete sim1;
    delete den2;
    delete den1;
    delete et;
    MPI_Finalize();
    
    return 0;
    
}