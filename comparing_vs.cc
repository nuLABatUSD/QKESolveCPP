#include <cmath>
#include <iostream>
#include <ostream>
#include <fstream>

#include "constants.hh"
#include "arrays.hh"
#include "QKE_methods.hh"

using namespace std;

int main(int argc, char* argv[]){
    
    const std::string& times_output_file = std::string(argv[4]);
    const std::string& v_m_mag_output_file = std::string(argv[3]);
    const std::string& vz_output_file = std::string(argv[2]);
    const std::string& input_file = std::string(argv[1]);
    
    linspace_and_gl* epsilon = new linspace_and_gl(0,20,201,5);
    double* dens_vals = new double[epsilon->get_len()*8+2]();
    
    std::ifstream densfile;
    densfile.open(input_file);

    std::ofstream v_m_z;
    v_m_z.open(vz_output_file);

    std::ofstream v_m_mag;
    v_m_mag.open(v_m_mag_output_file);
    
    std::ofstream times;
    times.open(times_output_file);
    

    if (!densfile.is_open()) {
        std::cout << "Error opening density input file" << std::endl;
    }
    if (!v_m_z.is_open()) {
        std::cout << "Error opening v_m_z output file" << std::endl;
    }
    if (!v_m_mag.is_open()) {
        std::cout << "Error opening v_m_mag output file" << std::endl;
    }
    if(!times.is_open()) {
        std::cout << "Error opening times output file" << std::endl;
    }
    
    three_vector_for_QKE* dummy_v_vacuum = new three_vector_for_QKE();
    dummy_v_vacuum->v_vacuum(0,1,0);
    std::cout << "v_vacuum magnitude=" << dummy_v_vacuum->magnitude();
    three_vector_for_QKE* dummy_v_density = new three_vector_for_QKE();
    
    
    std::string line;
    while(std::getline(densfile, line)){
        std::string densval;
        std::string delimiter = ", ";

        size_t pos = 0;
        int i=0;
        while((pos = line.find(delimiter)) != std::string::npos){
            densval = line.substr(0, pos);
            if(i==0){
                times << std::stod(densval) << ", ";
            }
            //this takes care of first two elements being initial place and initial step
            if(i>1){
                dens_vals[i-2] = std::stod(densval);
            }
            line.erase(0, pos + delimiter.length());
            i++;
        }
        dens_vals[i-2] = std::stod(line);
        density* dens = new density(epsilon->get_len(), epsilon, dens_vals);
        
        dummy_v_density->v_density(epsilon, dens);
        
        v_m_z << dummy_v_density->get_value(2) << ", ";
        v_m_mag << dummy_v_density->magnitude() << ", ";
        
        
        delete dens;

    }
    densfile.close();
    v_m_z.close();
    v_m_mag.close();
    
    
    delete[] dens_vals;
    delete epsilon;
    
    
    return 0; 
}

