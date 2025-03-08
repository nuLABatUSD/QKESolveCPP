#include <cmath>
#include <iostream>
#include <ostream>
#include <fstream>
#include <chrono>

#include "constants.hh"
#include "arrays.hh"
#include "QKE_methods.hh"
#include "mpi.h"
#include "collision_MPI.hh"

using namespace std;

double df_e_dt_plus_df_mu_dt(linspace_and_gl*, density*, int, bool);
double n(linspace_and_gl*, density*, bool);
double dn_dt(linspace_and_gl*, density*, bool, double*);
double p(linspace_and_gl*, density*, bool);
double dp_dt(linspace_and_gl*, density*, bool, double*);
double ds_dt_over_s(linspace_and_gl*, density*, bool, double*);
double summed_df_dt(linspace_and_gl*, density*, int, bool, bool);

int main(int argc, char* argv[]){
    int myid, numprocs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    
    //const std::string& entropy_output_file = std::string(argv[4]);
    const std::string& numberdens_output_file = std::string(argv[3]);
    const std::string& energy_output_file = std::string(argv[2]);
    const std::string& input_file = std::string(argv[1]);
    
    linspace_and_gl* epsilon = new linspace_and_gl(0,20,201,5);
    double* dens_vals = new double[epsilon->get_len()*8+2]();
    
    std::ifstream densfile;
    densfile.open(input_file);

    std::ofstream numberdensityfile;
    numberdensityfile.open(numberdens_output_file);

    std::ofstream energyfile;
    energyfile.open(energy_output_file);
    /*
    std::ofstream entropyfile;
    entropyfile.open(entropy_output_file);
    */

    

    if (!densfile.is_open()) {
        std::cout << "Error opening density input file" << std::endl;
    }
    if (!numberdensityfile.is_open()) {
        std::cout << "Error opening number density output file" << std::endl;
    }
    if (!energyfile.is_open()) {
        std::cout << "Error opening energy density output file" << std::endl;
    }
    /*
    if (!entropyfile.is_open()) {
        std::cout << "Error opening entropy output file" << std::endl;
    }*/
    
    collision_MPI* col = new collision_MPI(myid, numprocs, epsilon);
    
    double num_dens;
    double energy_dens;
    int j = 0;
    std::string line;
    while(std::getline(densfile, line)){
        auto start = std::chrono::high_resolution_clock::now();
        std::string densval;
        std::string delimiter = ", ";

        size_t pos = 0;
        int i=0;
        while((pos = line.find(delimiter)) != std::string::npos){
            densval = line.substr(0, pos);
            //this takes care of first two elements being initial place and initial step
            if(i>1){
                dens_vals[i-2] = std::stod(densval);
            }
            line.erase(0, pos + delimiter.length());
            i++;
        }
        dens_vals[i-2] = std::stod(line);
        density* dens = new density(epsilon->get_len(), epsilon, dens_vals);
        
        num_dens = col->num_dens_sum_rule(dens);
        energy_dens = col->energy_dens_sum_rule(dens);
        
        if(myid==0){
            numberdensityfile << num_dens << ", ";
            energyfile << energy_dens << ", ";
            std::cout << "i am adding " << num_dens << " to the number dens file and " << energy_dens << " to the energy file" << std::endl;
        }
        
        
        delete dens;

        auto stop = std::chrono::high_resolution_clock::now();

        auto duration = duration_cast<std::chrono::milliseconds>(stop - start);
        double time_elapsed = duration.count()/1000.;
        if(myid==0){std::cout << "finished line " << j << " in " << time_elapsed << "seconds" << std::endl;}
        j++;

    }
    densfile.close();
    numberdensityfile.close();
    energyfile.close();
    //entropyfile.close();
    
    
    delete[] dens_vals;
    delete epsilon;
    
    
    return 0; 
}
