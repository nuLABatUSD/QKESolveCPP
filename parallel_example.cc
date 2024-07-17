#include <iostream>
#include "mpi.h"
#include "arrays.hh"
#include <cmath>

using std::cout;
using std::endl;

double f(double, double);


int main(int argc, char *argv[]){
    int myid, numprocs;
    int sender, anstype;
    MPI_Status* status;
    int numpoints;
    double myans1, myans2, myans3;
    
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    
    numpoints = 1000;
    linspace_and_gl* xvals = new linspace_and_gl(0,10,numpoints,0);
    double* yvals  = new double[numpoints];
    
    if(myid == 0){
        for(int i=0; i<numpoints; i++){
            yvals[i] = exp(-1 * xvals->get_value(i));
        }
    }
    
    MPI_Bcast(yvals, numpoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
   
    if(myid != 0){
        dep_vars* y_vals = new dep_vars(numpoints);
        for(int i=0; i<numpoints; i++){
            y_vals->set_value(i, yvals[i]);
        }
        myans1 = xvals->integrate(y_vals);
        MPI_Send(&myans1, 1, MPI_DOUBLE, 0, myid-1, MPI_COMM_WORLD);
        
        delete y_vals;
    }
    
    if(myid == 0){
        
        //cout << "i am processor " << myid << " and i got here" << endl;
        dep_vars* integral_one_results = new dep_vars(3);
        for(int i=0; i<3; i++){
            MPI_Recv(&myans1, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status);
            sender = status->MPI_SOURCE;
            anstype = status->MPI_TAG;
            integral_one_results->set_value(anstype, myans1);
        }
        cout << "integral one: " << endl;
        integral_one_results->print_all();
        delete integral_one_results;
        
        for(int i=0; i<numpoints; i++){
            yvals[i] = exp(-2 * xvals->get_value(i));
        }
    }
    MPI_Bcast(yvals, numpoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if(myid != 0){
        dep_vars* y_vals = new dep_vars(numpoints);
        for(int i=0; i<numpoints; i++){
            y_vals->set_value(i, yvals[i]);
        }
        myans2 = xvals->integrate(y_vals);
        MPI_Send(&myans2, 1, MPI_DOUBLE, 0, myid-1, MPI_COMM_WORLD);
        
        delete y_vals;
    }
    
    if(myid == 0){
        
        //cout << "i am processor " << myid << " and i got here" << endl;
        dep_vars* integral_two_results = new dep_vars(3);
        for(int i=0; i<3; i++){
            MPI_Recv(&myans2, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status);
            sender = status->MPI_SOURCE;
            anstype = status->MPI_TAG;
            integral_two_results->set_value(anstype, myans2);
        }
        cout << "integral two: " << endl;
        integral_two_results->print_all();
        delete integral_two_results;
        
        for(int i=0; i<numpoints; i++){
            yvals[i] = exp(-3 * xvals->get_value(i));
        }
    }
    
    MPI_Bcast(yvals, numpoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(myid != 0){
        dep_vars* y_vals = new dep_vars(numpoints);
        for(int i=0; i<numpoints; i++){
            y_vals->set_value(i, yvals[i]);
        }
        myans3 = xvals->integrate(y_vals);
        MPI_Send(&myans3, 1, MPI_DOUBLE, 0, myid-1, MPI_COMM_WORLD);
        
        delete y_vals;
    }
    
    if(myid == 0){
        
        //cout << "i am processor " << myid << " and i got here" << endl;
        dep_vars* integral_three_results = new dep_vars(3);
        for(int i=0; i<3; i++){
            MPI_Recv(&myans3, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status);
            sender = status->MPI_SOURCE;
            anstype = status->MPI_TAG;
            integral_three_results->set_value(anstype, myans3);
        }
        cout << "integral three: " << endl;
        integral_three_results->print_all();
        delete integral_three_results;
    
    }
    delete[] yvals;
    delete xvals;
    MPI_Finalize();
    
    
    
}