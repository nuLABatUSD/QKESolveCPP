#include <iostream>
#include "arrays.hh"
#include "QKE_methods.hh"
#include <ostream>
#include <fstream>

using std::cout;
using std::endl;
using namespace std;


/*
TO RUN:
g++ initialize_dens.cc array_methods.cc matrices.cc thermodynamics.cc QKE_methods.cc -std=c++11 -o makedens
./makedens xmin, xmax, numlin, numgl, T, A, B, output file name

*/

//inputs are xmin, xmax, numlin, numgl, T, A, B, output file name
int main(int argc, char *argv[])
{
    
    double xmin = std::atof(argv[1]);
    double xmax = std::atof(argv[2]);
    int numlin = std::atoi(argv[3]);
    int numgl = std::atoi(argv[4]);
    double T = std::atof(argv[5]);
    int A = std::atoi(argv[6]);
    int B = std::atoi(argv[7]);
    const std::string& output_file = std::string(argv[8]);
    
    
    linspace_and_gl* et = new linspace_and_gl(xmin, xmax, numlin, numgl);
    density* den1 = new density(et, A, B);
    den1->set_T(T);
    
    std::ofstream myfile;
    myfile.open(output_file);
    for(int j=0; j<den1->length()-1; j++){
        myfile << den1->get_value(j) << ", ";
    }
    myfile << den1->get_value(den1->length()-1) << endl;
    myfile.close();
    
    delete den1;
    delete et;
    return 0;
    
}

