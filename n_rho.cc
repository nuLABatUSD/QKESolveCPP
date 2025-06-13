#include <iostream>
#include <fstream>
#include <chrono>

#include "collision_MPI.hh"
#include "arrays.hh"
#include "QKE_methods.hh"
#include "mpi.h"

using std::cout;
using std::endl;
using std::ofstream;

using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
using std::chrono::duration_cast;

int main(int argc, char* argv[])
{
    int myid, numprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    cout << myid << endl;

    linspace_and_gl* eps = new linspace_and_gl(0., 20., 201, 5);
    density* thermal = new density(eps, 1, 8);
    thermal->set_T(16.);


    collision_MPI* coll = new collision_MPI(myid, numprocs, eps);
    
    dep_vars** results = new dep_vars*[4];
    dep_vars** results_z = new dep_vars*[4];
    for(int j = 0; j < 4; j++){
        results[j] = new dep_vars(eps->get_len());
        results_z[j] = new dep_vars(eps->get_len());
    }
    coll->collision_term(thermal, results, results_z);
    
    if(myid == 0){
        ofstream file("R_results.csv");
        ofstream file_z("Rz_results.csv");
	ofstream de("R_dens_results.csv");

	thermal->print_csv(de);
        
        auto start = high_resolution_clock::now();
        
        for(int i = 0; i < eps->get_len(); i++){
            file << eps->get_value(i) << ", " << results[0]->get_value(i) << ", " << results[2]->get_value(i) << ", " << results[0]->get_value(i) / results[2]->get_value(i) << ", ";
            file << results[1]->get_value(i) << ", " << results[3]->get_value(i) << ", " << results[1]->get_value(i) / results[3]->get_value(i) << endl;
            
            file_z << eps->get_value(i) << ", " << results_z[0]->get_value(i) << ", " << results_z[2]->get_value(i) << ", " << results_z[0]->get_value(i) / results_z[2]->get_value(i) << ", ";
            file_z << results_z[1]->get_value(i) << ", " << results_z[3]->get_value(i) << ", " << results_z[1]->get_value(i) / results_z[3]->get_value(i) << endl;
            
        }
        
        file.close();
	file_z.close();
        de.close();

        dep_vars** ddt = new dep_vars*[4];
        for(int j = 0; j < 4; j++)
            ddt[j] = new dep_vars(eps->get_len());
            
        for(int i = 0; i < eps->get_len(); i++)
        {
            ddt[0]->set_value(i, pow(eps->get_value(i),2) * results[0]->get_value(i));
            ddt[1]->set_value(i, pow(eps->get_value(i),2) * results[2]->get_value(i));
            
            ddt[2]->set_value(i, pow(eps->get_value(i),3) * results[1]->get_value(i));
            ddt[3]->set_value(i, pow(eps->get_value(i),3) * results[3]->get_value(i));
        }
        
        cout << "Sum rule dn/dt = " << eps->integrate(ddt[0]) / eps->integrate(ddt[1]) << endl;
        cout << "Sum rule drho/dt = " << eps->integrate(ddt[2]) / eps->integrate(ddt[3]) << endl;
        
        for(int j = 0; j < 4; j++)
            delete ddt[j];
        delete[] ddt;
        
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        
        cout << duration.count() / 1000. << " seconds" << endl;
    }
    for(int j = 0; j < 4; j++)
        delete results[j];
    delete[] results;
    delete coll;
    delete eps;

    cout << "The end: " << myid << endl;

	MPI_Finalize();
    return 0;
}
