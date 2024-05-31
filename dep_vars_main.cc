#include <iostream>
#include <cmath>
#include "arrays.hh"
#include "constants.hh" // so now I can just write the names of variables in constants.hh


using std::cout;
using std::endl;

void f(double, density*, density*);

/*
    void f(double x, double* y, int N, double* der) 
{
    if (N == 3)
    {
        der[0] = 1;
        der[1] = x;
        der[2] = 1/(x+1);

    } else if (N == 4) {
        der[0] = y[1];
        der[1] = y[2];
        der[2] = y[3];
        der[3] = pow(_PI_, 4) * y[0]; // really its std::pow?
    } 
} 


void f(double x, dep_vars* y, dep_vars* z)
{
    int N = y->length();

    if (N==3)
    {
        z->set_value(0, 1.);
        z->set_value(1, x);
        z->set_value(2, 1/(1+x));
    }
    else if (N==4)
    {
        z->set_value(0, y->get_value(1));
        z->set_value(1, y->get_value(2));
        z->set_value(2, y->get_value(3));
        z->set_value(3, pow(_PI_, 4) * y->get_value(0));
    }
    return;
}

*/

int main()
{
   
   
    linspace_for_trap* et = new linspace_for_trap(0.,20, 201);
    three_vector* v = new three_vector();
    
    double eta_e = 0.2;
    double eta_mu = -0.02;
    

    density* den = new density(et, eta_e, eta_mu);
    
    density* new_den = new density(et->get_len(), et);

    f(0, den, new_den);
    new_den->print_all();
    
    
    delete et;
    delete den;
    delete v;
    return 0;
}

void f(double t, density* d1, density* d2){
    dummy_vars* E = d1->get_E();
    three_vector* dummy_v_vac = new three_vector();
    three_vector* dummy_v_dens = new three_vector();
    three_vector* dummy_v_therm = new three_vector();
    three_vector* dummy_v_1 = new three_vector();
    three_vector* dummy_v_2 = new three_vector();
    three_vector* dummy_v_3 = new three_vector();
    three_vector* vcrossp = new three_vector();
    
    three_vector* p = new three_vector();
    
    dummy_v_vac->v_vacuum();
    dummy_v_dens->v_density(E, d1);
    dummy_v_therm->v_thermal(E, d1);
    
    
    for (int i=1; i< E->get_len(); i++){
        dummy_v_1->copy(dummy_v_vac);
        dummy_v_1->multiply_by(1./E->get_value(i));

        dummy_v_2->copy(dummy_v_therm);
        dummy_v_2->multiply_by(E->get_value(i));
        
        dummy_v_3->add(dummy_v_vac, dummy_v_dens);
        dummy_v_3->add(dummy_v_1, dummy_v_therm);       

        
        d1->p_vector(i,true,p);
        
        vcrossp->set_cross_product(dummy_v_3,p);
        
        
        d2->set_value(4*i+1, vcrossp->get_value(0));
        d2->set_value(4*i+2, vcrossp->get_value(1));
        d2->set_value(4*i+3, vcrossp->get_value(2));
        
        
        d1->p_vector(i,false,p);
        vcrossp->set_cross_product(dummy_v_3,p);
        
        
        d2->set_value(4*(E->get_len())+4*i+1, vcrossp->get_value(0));
        d2->set_value(4*(E->get_len())+4*i+2, vcrossp->get_value(1));
        d2->set_value(4*(E->get_len())+4*i+3, vcrossp->get_value(2));
        
    }
        
    delete dummy_v_vac;
    delete dummy_v_dens;
    delete dummy_v_therm;
    delete dummy_v_1;
    delete dummy_v_2;
    delete dummy_v_3;
    delete vcrossp;
    delete p;
    
}