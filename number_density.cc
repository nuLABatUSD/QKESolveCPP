#include <cmath>
#include <iostream>
#include <ostream>
#include <fstream>
#include "constants.hh"

#include "arrays.hh"
#include "QKE_methods.hh"

using namespace std;

double df_e_dt_plus_df_mu_dt(linspace_and_gl*, density*, int, bool);
double n(linspace_and_gl*, density*, bool);
double dn_dt(linspace_and_gl*, density*, bool);
double p(linspace_and_gl*, density*, bool);
double dp_dt(linspace_and_gl*, density*, bool);

int main(int argc, char* argv[]){
    const std::string& numberdens_output_file = std::string(argv[3]);
    const std::string& energy_output_file = std::string(argv[2]);
    const std::string& input_file = std::string(argv[1]);
    
    linspace_and_gl* epsilon = new linspace_and_gl(0,10,201,5);
    double* dens_vals = new double[epsilon->get_len()*8+2]();

    
    std::ofstream numberdensityfile;
    numberdensityfile.open(energy_output_file);
    
    std::ofstream energyfile;
    energyfile.open(energy_output_file);
    
    std::ifstream densfile;
    densfile.open(input_file);
    
    
    if (!densfile.is_open()) {
        std::cout << "Error opening density input file" << std::endl;
    }
    if (!numberdensityfile.is_open()) {
        std::cout << "Error opening output file" << std::endl;
    }

    int j = 0;
    std::string line;
    while(std::getline(densfile, line)){
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
        double dn_dt_n = dn_dt(epsilon, dens, true) / n(epsilon, dens, true);
        double dp_dt_p = dp_dt(epsilon, dens, true) / p(epsilon, dens, true);
        numberdensityfile << dn_dt_n / (pow(dens->get_T(),5) * pow(_GF_,2)) << ", ";
        energyfile << dp_dt_p / (pow(dens->get_T(),5) * pow(_GF_,2)) << ", ";
        std::cout << "finished line " << j << std::endl;
        j++;
        
    }
    densfile.close();
    numberdensityfile.close();
    energyfile.close();
    
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

double dn_dt(linspace_and_gl* eps, density* dens, bool neutrino){
    dep_vars* int_vals = new dep_vars(eps->get_len());
    
    for(int i=0; i<eps->get_len(); i++){
        int_vals->set_value(i, pow(eps->get_value(i),2) * (df_e_dt_plus_df_mu_dt(eps, dens, i, neutrino)));
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

double dp_dt(linspace_and_gl* eps, density* dens, bool neutrino){
    dep_vars* int_vals = new dep_vars(eps->get_len());
    
    for(int i=0; i<eps->get_len(); i++){
        int_vals->set_value(i, pow(eps->get_value(i),3) * (df_e_dt_plus_df_mu_dt(eps, dens, i, neutrino)));
    }
    
    double result = eps->integrate(int_vals);
    delete int_vals;
    return result;
}

double df_e_dt_plus_df_mu_dt(linspace_and_gl* eps, density* dens, int i, bool neutrino){
    nu_nu_collision* nu_nu = new nu_nu_collision(eps, i);
    double* nu_nu_int = new double[4]();
    nu_nu->whole_integral(dens, neutrino, nu_nu_int);
    
    double Tcm = dens->get_T();
    nu_e_collision* nu_e = new nu_e_collision(eps, i, Tcm);
    double* nu_e_int = new double[4]();
    nu_e->whole_integral(dens, neutrino, nu_e_int);
    
    double result = nu_e_int[0] + nu_nu_int[0];
    delete nu_e;
    delete nu_nu;
    delete[] nu_nu_int;
    delete[] nu_e_int;
    return result;
    
}

