#include "arrays.hh"
#include <cmath>
#include <iostream>


using std::cout;
using std::endl;
using std::ostream;

//dummy_vars
dummy_vars::dummy_vars(int num){
    N = num;
    values = new double[N];
    weights = new double[N-1];
}

void dummy_vars::print_all(){
    for(int i =0; i<N; i++){
        cout << values[i] << endl;
    }
}

double dummy_vars::get_value(int i){
    return values[i];
}

void dummy_vars::set_value(int i, double v)
{values[i] = v;}

double dummy_vars::get_dx_val(int i){
    return weights[i];
}

int dummy_vars::get_len(){
    return N;
}

double dummy_vars::integrate(dep_vars* fvals){
    double result = 0;
    for (int i = 0; i<N; i++){
       result += fvals->get_value(i) * weights[i]; 
    }
    return result;    
}

dummy_vars::~dummy_vars(){
    delete[] values;
    delete[] weights;
}


//linspace
linspace::linspace(double xmin, double xmax, int num):dummy_vars(num)
{
    double dx_val = (xmax - xmin) / (N-1);
    for (int i = 0; i<N; i++){
        values[i] = xmin + dx_val * i;
    }
        
}

//linspace_for_trap
linspace_for_trap::linspace_for_trap(double xmin, double xmax, int num):linspace(xmin, xmax, num)
{
    double dx_val = (xmax - xmin) / (N-1);
    weights[0] = dx_val / 2;
    weights[N-1] = dx_val / 2;
    for (int i=1; i<N-2; i++){
        weights[i] = dx_val;
    }
}
