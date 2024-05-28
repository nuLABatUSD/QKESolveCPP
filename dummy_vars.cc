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
    dx = new double[N-1];
}

void dummy_vars::print_all(){
    for(int i =0; i<N; i++){
        cout << values[i] << endl;
    }
}

dummy_vars::~dummy_vars(){
    delete[] values;
    delete[] dx;
}


//linspace
linspace::linspace(double xmin, double xmax, int num):dummy_vars(num)
{

    double dx_val = (xmax - xmin) / (N-1);
    for (int i = 0; i<N; i++){
        dx[i] = dx_val;
        values[i] = xmin + dx_val * i;
    }
        
}
