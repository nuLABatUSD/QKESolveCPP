#include <iostream>

#include "QKE_methods.hh"
#include "arrays.hh"
#include "collision_MPI.hh"
#include "mpi.h"

using std::cout;
using std::endl;

collision_MPI::collision_MPI(int rank, int num_ranks, linspace_and_gl* e)
{
    myid = rank;
    numprocs = num_ranks;
    
    eps = new linspace_and_gl(e);
    
    int_objects = new nu_nu_collision*[eps->get_len()];
    if (myid != 0)
        for (int i = myid-1; i < eps->get_len(); i += numprocs-1)
            int_objects[i] = new nu_nu_collision(eps, i);
}

collision_MPI::~collision_MPI(){
    if(myid != 0)
        for(int i = myid-1; i < eps->get_len(); i += numprocs-1)
            delete int_objects[i];
    delete[] int_objects;
    
    delete eps;
}

void collision_MPI::calculate_R(density* input, dep_vars** output){
    if (input->num_bins() != output[0]->length()){
        cout << "collision_MPI::calculate_R the input and output objects need the same length" << endl;
        return;
    }
    
    double* out_vals = new double[4 * input->num_bins()];
    double my_ans = 0.;
    int sender, tag;
    double* dummy_int = new double[4];
    MPI_Status status;
    
    if(myid == 0){
        for(int i = 0; i < eps->get_len(); i++){
            MPI_Recv(dummy_int, 4, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            sender = status.MPI_SOURCE;
            tag = status.MPI_TAG;
            
            for(int j = 0; j < 4; j++)
                out_vals[4*tag+j] += dummy_int[j];
            
            cout << "Received: " << tag << endl;
        }
        
    }
    else{
        double** nu_nu_int = new double*[4];
        for (int j = 0; j < 4; j++)
            nu_nu_int[j] = new double[4]();
        
        for(int i = myid-1; i < eps->get_len(); i += numprocs-1){
            int_objects[i]->whole_integral(input, true, nu_nu_int[0], true);
            int_objects[i]->whole_integral(input, false, nu_nu_int[1], true);
            
            int_objects[i]->whole_integral(input, true, nu_nu_int[2], false);
            int_objects[i]->whole_integral(input, false, nu_nu_int[3], false);
            
            for(int j = 0; j < 4; j++)
                dummy_int[j] = nu_nu_int[j][0];
            MPI_Send(dummy_int, 4, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
        }
        
        for(int j = 0; j < 4; j++)
            delete nu_nu_int[j];
        delete[] nu_nu_int;
    }
    
    MPI_Bcast(out_vals, 4 * input->num_bins(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    for(int j = 0; j < 4; j++)
        for(int i = 0; i < input->num_bins(); i++)
            output[j]->set_value(i, out_vals[4*i+j]);
            
    delete[] out_vals;
    delete[] dummy_int;
    
}

void collision_MPI::collision_term(density* input, dep_vars** output_0, dep_vars** output_z){
    if (input->num_bins() != output_0[0]->length()){
        cout << "collision_MPI::collision_term the input and output objects need the same length" << endl;
        return;
    }
    
    double* out_vals = new double[8 * input->num_bins()];
    double my_ans = 0.;
    int sender, tag;
    double* dummy_int = new double[8];
    MPI_Status status;
    
    if(myid == 0){
        for(int i = 0; i < eps->get_len(); i++){
            MPI_Recv(dummy_int, 8, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            sender = status.MPI_SOURCE;
            tag = status.MPI_TAG;
            
            for(int j = 0; j < 8; j++){
                out_vals[8*tag+j] = dummy_int[j];
            }
            cout << "Received: " << tag << endl;
        }
        
    }
    else{
        double** nu_nu_int = new double*[4];
        for (int j = 0; j < 4; j++)
            nu_nu_int[j] = new double[4]();
        
        for(int i = myid-1; i < eps->get_len(); i += numprocs-1){
            int_objects[i]->whole_integral(input, true, nu_nu_int[0], true);
            int_objects[i]->whole_integral(input, false, nu_nu_int[1], true);
            
            int_objects[i]->whole_integral(input, true, nu_nu_int[2], false);
            int_objects[i]->whole_integral(input, false, nu_nu_int[3], false);
            
            for(int j = 0; j < 4; j++){
                dummy_int[j] = nu_nu_int[j][0];
                dummy_int[4+j] = nu_nu_int[j][3];
            }
            MPI_Send(dummy_int, 8, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
        }
        
        for(int j = 0; j < 4; j++)
            delete nu_nu_int[j];
        delete[] nu_nu_int;
    }
    
    MPI_Bcast(out_vals, 8 * input->num_bins(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    for(int j = 0; j < 4; j++)
        for(int i = 0; i < input->num_bins(); i++){
            output_0[j]->set_value(i, out_vals[8*i+j]);
            output_z[j]->set_value(i, out_vals[8*i+4+j]);
        }
            
    delete[] out_vals;
    delete[] dummy_int;

}
