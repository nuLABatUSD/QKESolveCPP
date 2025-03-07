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


//input parameters will be xmin (1), xmax (2), numlin (3), numgl (4), sin2theta (5), deltamsquared (6), N_step (7), dN (8), x_initial (9), x_final (10), dx_initial (11), verbose (12), output file name (13), density file name (14)
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
    
    double sin2theta = std::atof(argv[5]);
    double deltamsquared = std::atof(argv[6]);
    
    int N_step = std::atoi(argv[7]);
    int dN = std::atoi(argv[8]);
    double x_0 = std::atof(argv[9]);
    double x_f = std::atof(argv[10]);
    double dx_0 = std::atof(argv[11]);
    bool verbose = argv[12];
    
    
    const std::string& output_file = std::string(argv[13]);
    const std::string& input_file = std::string(argv[14]);
    
    
    linspace_and_gl* et = new linspace_and_gl(xmin, xmax, numlin, numgl);
    QKESolveMPI* sim1 = new QKESolveMPI(myid, numprocs, et, sin2theta, deltamsquared, x_0, dx_0, input_file);
    
    sim1->just_neutrino_collision();
    
    
    
//int N_step, int dN, double x_final, const std::string& file_name, bool verbose = false
    sim1->run(N_step, dN, x_f, output_file, verbose);

    delete sim1;
    delete et;
    MPI_Finalize();
    return 0;
    
}

