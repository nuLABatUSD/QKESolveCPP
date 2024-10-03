#include <cmath>
#include <iostream>
#include <ostream>
#include <fstream>
#include <chrono>

#include "constants.hh"
#include "arrays.hh"
#include "QKE_methods.hh"

using namespace std;

double df_e_dt_plus_df_mu_dt(linspace_and_gl*, density*, int, bool);
double n(linspace_and_gl*, density*, bool);
double dn_dt(linspace_and_gl*, density*, bool, double*);
double p(linspace_and_gl*, density*, bool);
double dp_dt(linspace_and_gl*, density*, bool, double*);
double ds_dt_over_s(linspace_and_gl*, density*, bool, double*);

int main(int argc, char* argv[]){
    const std::string& entropy_output_file = std::string(argv[4]);
    const std::string& numberdens_output_file = std::string(argv[3]);
    const std::string& energy_output_file = std::string(argv[2]);
    const std::string& input_file = std::string(argv[1]);
    
    linspace_and_gl* epsilon = new linspace_and_gl(0,10,201,5);
    double* dens_vals = new double[epsilon->get_len()*8+2]();

    std::ofstream numberdensityfile;
    numberdensityfile.open(numberdens_output_file);

    std::ofstream energyfile;
    energyfile.open(energy_output_file);

    std::ofstream entropyfile;
    entropyfile.open(entropy_output_file);

    std::ifstream densfile;
    densfile.open(input_file);


    if (!densfile.is_open()) {
        std::cout << "Error opening density input file" << std::endl;
    }
    if (!numberdensityfile.is_open()) {
        std::cout << "Error opening number density output file" << std::endl;
    }
    if (!energyfile.is_open()) {
        std::cout << "Error opening energy density output file" << std::endl;
    }
    if (!entropyfile.is_open()) {
        std::cout << "Error opening entropy output file" << std::endl;
    }

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
        bool neutrino = true;
        double* C0_vals = new double[epsilon->get_len()]();
        
        
        for(int k=0; k<epsilon->get_len(); k++){
            C0_vals[k] = df_e_dt_plus_df_mu_dt(epsilon, dens, k, neutrino);
        }
        
        
        double multiplicative_factor = 1 / (pow(dens->get_T(),5) * pow(_GF_,2));
        std::cout << "multiplicative factor=" << multiplicative_factor << std::endl;
        double dn_dt_n = dn_dt(epsilon, dens, neutrino, C0_vals) / n(epsilon, dens, neutrino) * multiplicative_factor;
        //double dp_dt_p=0;
        double ds_dt_s=0;
        double dp_dt_p = dp_dt(epsilon, dens, neutrino, C0_vals) / p(epsilon, dens, neutrino) * multiplicative_factor;
        //double ds_dt_s = ds_dt_over_s(epsilon, dens, neutrino, C0_vals) * multiplicative_factor;

        std::cout << n(epsilon, dens, neutrino) << std::endl;
        numberdensityfile << std::to_string(dn_dt_n) << ", ";
        energyfile << std::to_string(dp_dt_p) << ", ";
        entropyfile << std::to_string(ds_dt_s) << ", ";
        std::cout << "i am adding " << dn_dt_n << " to the number dens file, " << dp_dt_p << " to the energy file, and " << ds_dt_s << " to the entropy file" << std::endl;
        
        delete[] C0_vals;
        delete dens;

        auto stop = std::chrono::high_resolution_clock::now();

        auto duration = duration_cast<std::chrono::milliseconds>(stop - start);
        double time_elapsed = duration.count()/1000.;
        std::cout << "finished line " << j << " in " << time_elapsed << "seconds" << std::endl;
        j++;

    }
    densfile.close();
    numberdensityfile.close();
    energyfile.close();
    entropyfile.close();
    
    
    delete[] dens_vals;
    delete epsilon;
    
    
    return 0; 
}

double n(linspace_and_gl* eps, density* dens, bool neutrino){
    dep_vars* int_vals = new dep_vars(eps->get_len());
    
    for(int i=0; i<eps->get_len(); i++){
        int_vals->set_value(i, pow(eps->get_value(i),2) * dens->p0(i, neutrino));
    }
    
    double result = eps->integrate(int_vals);
    delete int_vals;
    return result;
}

double dn_dt(linspace_and_gl* eps, density* dens, bool neutrino, double* C0vals){
    dep_vars* int_vals = new dep_vars(eps->get_len());
    
    for(int i=0; i<eps->get_len(); i++){
        int_vals->set_value(i, pow(eps->get_value(i),2) * C0vals[i]);
    }
    double result = eps->integrate(int_vals);
    delete int_vals;
    return result;
}

double p(linspace_and_gl* eps, density* dens, bool neutrino){
    dep_vars* int_vals = new dep_vars(eps->get_len());
    
    for(int i=0; i<eps->get_len(); i++){
        int_vals->set_value(i, pow(eps->get_value(i),3) * dens->p0(i, neutrino));
    }
    
    double result = eps->integrate(int_vals);
    delete int_vals;
    return result;
}

double dp_dt(linspace_and_gl* eps, density* dens, bool neutrino, double* C0vals){
    dep_vars* int_vals = new dep_vars(eps->get_len());
    
    for(int i=0; i<eps->get_len(); i++){
        int_vals->set_value(i, pow(eps->get_value(i),3) * C0vals[i]);
    }
    
    double result = eps->integrate(int_vals);
    delete int_vals;
    return result;
}

double ds_dt_over_s(linspace_and_gl* eps, density* dens, bool neutrino, double* C0vals){
    dep_vars* s_int_vals = new dep_vars(eps->get_len());
    dep_vars* ds_dt_int_vals = new dep_vars(eps->get_len());
    
    for(int i=0; i<eps->get_len(); i++){
        double logp0 = log(dens->p0(i, neutrino));
        s_int_vals->set_value(i, pow(eps->get_value(i),2) * logp0);
        ds_dt_int_vals->set_value(i, pow(eps->get_value(i),2) * (logp0 + 1) * C0vals[i]);
    }
    
    double s = eps->integrate(s_int_vals);
    double ds_dt = eps->integrate(ds_dt_int_vals);
    
    delete s_int_vals;
    delete ds_dt_int_vals;
    
    return ds_dt / s;
}

double df_e_dt_plus_df_mu_dt(linspace_and_gl* eps, density* dens, int i, bool neutrino){
    nu_nu_collision* nu_nu = new nu_nu_collision(eps, i);
    double* nu_nu_int = new double[4]();
    nu_nu->whole_integral(dens, neutrino, nu_nu_int);
    
    double Tcm = dens->get_T();
    //std::cout << "eps->length is " << eps->get_len() << " and i am using index " << i << std::endl;
    nu_e_collision* nu_e = new nu_e_collision(eps, i, Tcm);
    double* nu_e_int = new double[4]();
    nu_e->whole_integral(dens, neutrino, nu_e_int);
    
    //std::cout << "nu_nu=" << nu_nu_int[0] << ", nu_e=" << nu_e_int[0] << ", ";
    double result = nu_e_int[0] + nu_nu_int[0];
    delete nu_e;
    delete nu_nu;
    delete[] nu_nu_int;
    delete[] nu_e_int;
    return result;
    
}

