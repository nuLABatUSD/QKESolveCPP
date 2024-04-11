#include <iostream>
#include "quad_vals.hh"
#include <cmath>

double f(double u, double m, double T);
double integral(double m, double T);

using std::cout;

int main(){
    double T = 100;
    double m = 0.511;
    
    cout << integral(m,T);
}

double f(double u, double m, double T){
    return pow(u,2) * pow(T,3) * sqrt(pow(u,2)*pow(T,2) + pow(m,2)) / (exp(sqrt(pow(u,2)*pow(T,2) + pow(m,2))/T) + 1) * exp(u);
}                                                             
                                                               
double integral(double m, double T){
    double sum = 0;
    for (int i=0; i<nvals; i++){
        sum += f(xvals[i],m,T) * wvals[i];
    }
    return sum;
}