#include <iostream>
#include "QKESolveMPI.hh"
#include "arrays.hh"
#include "QKE_methods.hh"
#include <chrono>
#include <fstream>


#include "mpi.h"
using std::cout;
using std::endl;
using namespace std;


/*
TO RUN:
mpic++ test4.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -std=c++11 -o wed
mpiexec -n 4 wed same "outputs.csv" "times.csv"

*/

int main(int argc, char *argv[])
{
    int myid, numprocs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    linspace_and_gl* et = new linspace_and_gl(0., 10., 201, 5);
    double eta_e = 0.01;
    double eta_mu = -0.01;
    std::ofstream resultsfile;
    if(myid==0){
        
        resultsfile.open(argv[2]);}
    
    
    double* time_vals = new double[21]();

    for(int i=-10; i<11; i++){
        eta_e = i * 0.01;
        
        if(strcmp(argv[1],"same")==0){
            eta_mu = eta_e;
        }
        else if(strcmp(argv[1],"opposite")==0){
            eta_mu = -1 * eta_e;
        }
        else if(strcmp(argv[1],"double")==0){
            eta_mu = 2 * eta_e;
        }
        else if(strcmp(argv[1],"half")==0){
            eta_mu = .5 * eta_e;
        }
        else{
            cout << "this type not supported!" << endl;
        }
        
        QKESolveMPI* sim1 = new QKESolveMPI(myid, numprocs, et, 0.8, 2.5e-15, eta_e, eta_mu);
        density* den1 = new density(et, eta_e, eta_mu);
        density* den2 = new density(den1);
        den1->set_T(0.25);
        sim1->set_ics(0, den1, 1.e12);



        auto start = high_resolution_clock::now();
        sim1->f(1, den1, den2);

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        double time_elapsed = duration.count()/1000.;

        double max_time_elapsed = 0;
        MPI_Reduce(&time_elapsed, &max_time_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        
        
        if(myid == 0){
            //prints to csv where first value is eta_e, second is eta_mu
            //each row in csv contains one whole derivative
            
            resultsfile << eta_e << ", " << eta_mu << ", ";
            for(int j=0; j<den2->length()-1; j++){
                resultsfile << den2->get_value(j) << ", ";
            }
            resultsfile << den2->get_value(den2->length()-1) << endl;
            
            
            
        }
        time_vals[i] = max_time_elapsed;
        

        
        delete sim1;
        delete den1;
        delete den2;
        
        
    }
    
    if(myid==0){
        resultsfile.close();

        std::ofstream timefile;
        timefile.open(argv[3]);
        for(int j=0; j<21; j++){
            timefile << time_vals[j] << endl;
        }
    timefile.close();}
    delete[] time_vals;
    delete et;
    MPI_Finalize();
    return 0;
    
}

