#include <iostream>
#include <chrono>

using std::cout;
using std::endl;
using namespace std::chrono;

#include "arrays.hh"
#include "QKE_methods.hh"

// g++ profile_integrals.cc array_methods.cc QKE_methods.cc matrices.cc thermodynamics.cc -std=c++11 -o integ
// sample ##### 30 -file sample.txt  << enter process number for ######
// run on up_quark: elapsed time (single processor) 379 sec (6m:19s)

int main()
{
    linspace_and_gl* et = new linspace_and_gl(0., 10., 201, 5);
    int N_eps = et->get_len();
    double eta_e = 0.01;
    double eta_mu = -0.01;

    density* den1 = new density(et, eta_e, eta_mu);
    den1->set_T(0.25);

    integration** integrators = new integration*[N_eps];
    for(int i = 0; i < N_eps; i++)
        integrators[i] = new integration(et, i);

    density* den2 = new density(den1);
    den2->zeros();

    auto start = high_resolution_clock::now();

    for(int i = 0; i < N_eps; i++)
        { cout << i << endl;
        for(int j = 0; j < 4; j++)
            {   
                den2->set_value(i*4 + j + 1, integrators[i]->whole_integral(den1, true, j));
                
                den2->set_value(N_eps*4 + i*4 + j + 1, integrators[i]->whole_integral(den1, false, j));
            }
        }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

        cout << endl << "Time elapsed: "
     << duration.count()/1000. << " seconds" << endl;


    delete den1;
    delete den2;

    for (int i = 0; i < N_eps; i++)
        delete integrators[i];
    delete[] integrators;
    
    return 0;
}