#include <iostream>
#include <cmath>
#include "arrays.hh"
#include "constants.hh" // so now I can just write the names of variables in constants.hh


using std::cout;
using std::endl;

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
*/

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

int main()
{
   
   
    linspace* et = new linspace(0.,20, 201);
    three_vector* v = new three_vector();
    
    double eta_e = 0.1;
    double eta_mu = -0.01;
    
    double v_dens = 0;
    v_dens += pow(_PI_,2) / 3 * (eta_e - eta_mu) + 1/3 * (pow(eta_e,3)-pow(eta_mu,3));

    
    double v_therm = 0;
    v_therm += pow(_PI_,2) / 2 * (pow(eta_e,2) - pow(eta_mu,2)) + 1/4 * (pow(eta_e,4) - pow(eta_mu,4));

    
    
    cout << "expected dens: " << v_dens << endl;
    cout << "expected therm: " << v_therm << endl;
    
    

    density* den = new density(et, eta_e, eta_mu);
    
    
    v->v_density(et, den);
    v->print_all();

   
    delete et;
    delete den;
    delete v;
    return 0;
}