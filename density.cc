#include <iostream>
#include "arrays.hh"

using std::cout;
using std::endl;

density::density(int num):dep_vars(num)
{
    N_bins = num;
    N = 8 * N_bins;
    dens = new dep_vars(N);
}

density::density(linspace eps, double eta_nu, double eta_mu):dep_vars(4)
{
    //N_bins = eps->length();
    N_bins = 5; 
    N = 8 * N_bins;
    dens = new dep_vars(N);
    //for (int i=0; i<N; i++){
     //   dens->set_value(i, 3);
    //}
    
}

int density::num_bins()
{return N_bins;}

void density::p_vector(int t, bool neutrino, three_vector* p)
{
    p = new three_vector();
    if(neutrino==true){
        for(int i=0; i<3; i++){
            p->set_value(i, dens->get_value(t+i+1));
        }}
    else{
        for(int i=0; i<3; i++){
            p->set_value(i, dens->get_value(N/2+t+i));
        }}
}

void density::p0_p(int t, bool neutrino, three_vector* p)
{
    if(neutrino==true){
        t=5;
    }
}