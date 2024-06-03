//#include "matrices.hh"
#include <iostream>
#include <complex> 
#include "matrices.hh"

using std::cout;
using std::endl;
using std::complex;

int main(){
    complex<double> z (3.2, 4.1);
    complex<double> x (1.0,1.0);
    
    complex_three_vector* p = new complex_three_vector(z);
    complex_three_vector* q = new complex_three_vector(x);
    complex_three_vector* s = new complex_three_vector(3);
    complex<double> y = q->magnitude_squared();
    cout << y << endl;
    
    delete p;
    delete q;
    delete s;
    
    return 0;
}