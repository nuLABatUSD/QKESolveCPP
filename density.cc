#include <iostream>
#include "arrays.hh"
#include <cmath>

using std::cout;
using std::endl;

density::density(int num):dep_vars(8*num)
{
    N_bins = num;
}

density::density(linspace* eps, double eta_nu, double eta_mu):dep_vars(8*eps->N)
{
    N_bins = eps->N;

    double fnu = 0;
    double fmu = 0;
    double fnubar = 0;
    double fmubar = 0;
    
    for (int i=0; i<N_bins; i++){
        fnu = 1 / (exp(eps->values[i] - eta_nu)+1);
        fmu = 1 / (exp(eps->values[i] - eta_mu)+1);
        values[4*i] = fnu + fmu;
        values[4*i+3] =  (fnu - fmu)/(fnu+fmu);
        //cout << i << ": " << dens->get_value(4*i) << endl;
       
       fnubar = 1 / (exp(eps->values[i] + eta_nu)+1);
       fmubar = 1 / (exp(eps->values[i] + eta_mu)+1);
       values[4*N_bins + 4*i] = fnubar - fmubar;
       values[4*N_bins + 4*i+3] = (fnubar - fmubar)/(fnubar+fmubar);
      
    }
    
    
}

int density::num_bins()
{return N_bins;}

void density::p_vector(int t, bool neutrino, three_vector* p)
{
    if(neutrino==true){
        for(int i=0; i<3; i++){
            p->set_value(i, values[4*t+i+1]);
        }}
    else{
        for(int i=0; i<3; i++){
            p->set_value(i, values[N/2+4*t+1+i]);
        }}
}

void density::p0_p(int t, bool neutrino, three_vector* p)
{
    if(neutrino==true){
        for(int i=0; i<3; i++){
            p->set_value(i,values[4*t+i+1]);
        }
        p->multiply_by(values[t]);
        
    }
    
    else{
        for(int i=0; i<3; i++){
            p->set_value(i, values[N/2+4*t+i+1]);
        }
        p->multiply_by(values[t]);
    }
}