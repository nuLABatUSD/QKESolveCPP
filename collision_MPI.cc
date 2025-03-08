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
    //output_0 is (FRS/neutrino combo, energy bins), so 4 by 206
    //output_0[0]->neutrino, net
    //output_0[1]->antineutrino, net
    //output_0[2]->neutrino, FRS
    //output_0[3]->antineutrino, FRS
    
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

double collision_MPI::dn_dt(int index, dep_vars** C0vals){
    dep_vars* int_vals = new dep_vars(eps->get_len());
    
    for(int i=0; i<eps->get_len(); i++){
        int_vals->set_value(i, pow(eps->get_value(i),2) * C0vals[index]->get_value(i));
    }
    double result = eps->integrate(int_vals);
    delete int_vals;
    return result;
}

double collision_MPI::dp_dt(int index, dep_vars** C0vals){
    dep_vars* int_vals = new dep_vars(eps->get_len());
    
    for(int i=0; i<eps->get_len(); i++){
        int_vals->set_value(i, pow(eps->get_value(i),3) * C0vals[index]->get_value(i));
    }
    
    double result = eps->integrate(int_vals);
    delete int_vals;
    return result;
}

double collision_MPI::num_dens_sum_rule(density* dens){
    dep_vars** C0_vals = new dep_vars*[4];
    dep_vars** Cz_vals = new dep_vars*[4];
    for(int i=0; i<4; i++){
        C0_vals[i] = new dep_vars(eps->get_len());
        Cz_vals[i] = new dep_vars(eps->get_len());
    }
    
    this->collision_term(dens, C0_vals, Cz_vals);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    double dn_dt_net = 0;
    double dn_dt_FRS = 0;
    
    dn_dt_net = this->dn_dt(0, C0_vals) + this->dn_dt(1, C0_vals);
    dn_dt_FRS = this->dn_dt(2, C0_vals) + this->dn_dt(3, C0_vals);

    
    return dn_dt_net / dn_dt_FRS;
    
}

double collision_MPI::energy_dens_sum_rule(density* dens){
    dep_vars** C0_vals = new dep_vars*[4];
    dep_vars** Cz_vals = new dep_vars*[4];
    for(int i=0; i<4; i++){
        C0_vals[i] = new dep_vars(eps->get_len());
        Cz_vals[i] = new dep_vars(eps->get_len());
    }
    
    this->collision_term(dens, C0_vals, Cz_vals);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    double dp_dt_net = 0;
    double dp_dt_FRS = 0;
    
    dp_dt_net = this->dp_dt(0, C0_vals) + this->dp_dt(1, C0_vals);
    dp_dt_FRS = this->dp_dt(2, C0_vals) + this->dp_dt(3, C0_vals);

    
    
    return dp_dt_net / dp_dt_FRS;
    
}
