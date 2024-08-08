#include <iostream>
#include "QKESolveMPI.hh"
#include "arrays.hh"
#include "QKE_methods.hh"
#include <chrono>


#include "mpi.h"
using std::cout;
using std::endl;
using namespace std;


/*
TO RUN:
mpic++ test2.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -std=c++11 -o wed
mpiexec -n 4 wed

*/


//input parameters will be xmin (1), xmax (2), numlin (3), numgl (4), eta_e (5), eta_mu (6), sin2theta (7), deltamsquared (8), N_step (9), dN (10), x_initial (11), x_final (12), dx_initial (13), verbose (14), output file name (15), density file name (16)
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
    
    int N_step = std::atoi(argv[9]);
    int dN = std::atoi(argv[10]);
    double x_0 = std::atof(argv[11]);
    double x_f = std::atof(argv[12]);
    double dx_0 = std::atof(argv[13]);
    bool verbose = argv[14];
    
    const std::string& output_file = std::string(argv[15]);
    const std::string& input_file = std::string(argv[16]);
    
    
    linspace_and_gl* et = new linspace_and_gl(xmin, xmax, numlin, numgl);

    QKESolveMPI* sim1 = new QKESolveMPI(myid, numprocs, et, sin2theta, deltamsquared, eta_e, eta_mu, x_0, dx_0, input_file);
    
//int N_step, int dN, double x_final, const std::string& file_name, bool verbose = false
    sim1->run(N_step, dN, x_f, output_file, verbose);

    delete sim1;
    delete et;
    MPI_Finalize();
    return 0;
    
}

