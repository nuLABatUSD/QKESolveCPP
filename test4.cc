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

void make_initial_dens(double, double, int, int, double, double, double, const std::string&);

/*
TO RUN:
mpic++ test4.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -std=c++11 -o wed
mpiexec -n 4 wed 0. 10. 201 5 0.01 -0.01 0.8 2.5e-15 0. 1.e12 same "outputs.csv" "times.csv"

*/

//input parameters will be xmin (1), xmax (2), numlin (3), numgl (4), eta_e (5), eta_mu (6), sin2theta (7), deltamsquared (8), x_initial (9), dx_initial (10), string that indicates choice of eta_e/eta_mu relationship (11), results output file (12), times output file (13)
int main(int argc, char *argv[])
{
    int myid, numprocs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    double xmin = std::atof(argv[1]);
    double xmax = std::atof(argv[2]);
    int numlin = std::atoi(argv[3]);
    int numgl = std::atoi(argv[4]);
    
    double eta_e = std::atof(argv[5]);
    double eta_mu = std::atof(argv[6]);
    double sin2theta = std::atof(argv[7]);
    double deltamsquared = std::atof(argv[8]);
    
    double x_0 = std::atof(argv[9]);
    double dx_0 = std::atof(argv[10]);
       
    const std::string& results_output = std::string(argv[12]);
    const std::string& times_output = std::string(argv[13]);
    
    
    linspace_and_gl* et = new linspace_and_gl(xmin, xmax, numlin, numgl);
    
    std::ofstream resultsfile;
    if(myid==0){
        resultsfile.open(results_output);
    }
    
    double* time_vals = new double[21]();
    double num_trials = 21;
    int count = 0;

    for(int i=-10; i<11; i++){
        eta_e = i * 0.01;
        
        if(strcmp(argv[11],"same")==0){
            eta_mu = eta_e;
        }
        else if(strcmp(argv[11],"opposite")==0){
            eta_mu = -1 * eta_e;
        }
        else if(strcmp(argv[11],"double")==0){
            eta_mu = 2 * eta_e;
        }
        else if(strcmp(argv[11],"half")==0){
            eta_mu = .5 * eta_e;
        }
        else{
            cout << "this type not supported!" << endl;
        }
        
        
        make_initial_dens(xmin, xmax, numlin, numgl, 0.25, eta_e, eta_mu, "initialthermaldistribution.csv");
        QKESolveMPI* sim = new QKESolveMPI(myid, numprocs, et, sin2theta, deltamsquared, eta_e, eta_mu, x_0, dx_0, "initialthermaldistribution.csv");
        
        density* den1 = new density(et, eta_e, eta_mu);
        density* den2 = new density(den1);
        den1->set_T(0.25);


        auto start = high_resolution_clock::now();
        sim->f(1, den1, den2);
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
        
        time_vals[count] = max_time_elapsed;
        count++;
        
        delete sim;
        delete den1;
        delete den2;
    }
    
    if(myid==0){
        resultsfile.close();

        std::ofstream timefile;
        timefile.open(times_output);
        for(int j=0; j<num_trials; j++){
            timefile << time_vals[j] << endl;
        }
        timefile.close();
    }
    
    delete[] time_vals;
    delete et;
    MPI_Finalize();
    return 0;
    
}


void make_initial_dens(double xmin, double xmax, int numlin, int numgl, double T, double eta_e, double eta_mu, const std::string& output_file){
    
    linspace_and_gl* et = new linspace_and_gl(xmin, xmax, numlin, numgl);
    density* den1 = new density(et, eta_e, eta_mu);
    den1->set_T(T);
    
    std::ofstream myfile;
    myfile.open(output_file);
    for(int j=0; j<den1->length()-1; j++){
        myfile << den1->get_value(j) << ", ";
    }
    myfile << den1->get_value(den1->length()-1) << endl;
    myfile.close();
    
    delete den1;
    delete et;
    
}

