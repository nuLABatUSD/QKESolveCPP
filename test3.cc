#include <iostream>
#include "arrays.hh"
#include "QKE_methods.hh"
#include <chrono>
#include "alternative_integrals.hh"
#include <iomanip>
#include <fstream>

using std::cout;
using std::endl;
double interpolate_one(double, int, double*, double*);
double extrapolate_exponential_one(double, double, double, double, double);
double extrapolate_linear_one(double, double, double, double, double);
double interpolate_two(double, double, double, double, double);

int main(int argc, char *argv[])
{
    const std::string& term1_output = std::string(argv[1]);
    const std::string& term2_output = std::string(argv[2]);
    
    
    
    linspace_and_gl* eps = new linspace_and_gl(0,10,201,5);
    
    density* dens = new density(eps, 0.01, -0.01);
    dens->set_T(1);
    double* nu_nu_int = new double[4]();
    bool neutrino = true;
    
    int p1 = 75;
    int p2 = 205;
    double F0 = 0;
    three_vector* F = new three_vector();
    nu_nu_collision_one* og_collision = new nu_nu_collision_one(eps, p1);
    
    
    
    std::ofstream term1;
    term1.open(term1_output);

    std::ofstream term2;
    term2.open(term2_output);
    
    double interpolated_p0 = 0;
    for(int p3=eps->get_len()-1; p3>=0; p3--){
        /*
        og_collision->Fvvsc_components_term_1(dens, neutrino, p2, p3, &F0, F);
        term1 << F0 << ", ";
        og_collision->Fvvsc_components_term_2(dens, neutrino, p2, p3, &F0, F);
        term2 << F0 << ", ";*/
        double p4_energy = eps->get_value(p1)+ eps->get_value(p2) - eps->get_value(p3);
        int count = eps->index_below_for_interpolation(p4_energy);
        if (p4_energy < eps->get_value(eps->get_len()-1)){
            interpolated_p0 = interpolate_two(p4_energy, eps->get_value(count), eps->get_value(count+1), dens->p0(count, neutrino), dens->p0(count+1, neutrino));
        }
        else{
            //std::cout << "extrapolating for p4_energy=" << p4_energy << std::endl;
            interpolated_p0 = extrapolate_exponential_one(p4_energy, eps->get_value(count-2), eps->get_value(count-1), dens->p0(count-2, neutrino), dens->p0(count-1, neutrino));
        }
        term1 << interpolated_p0 << ", ";
        term2 << std::setprecision(16) << 1-interpolated_p0 << ",";
        
    }
    
    
    term1.close();
    term2.close();
        
    delete F;
    delete og_collision;
    
    
    /*
    for(int p1=74; p1<75; p1++){
        nu_nu_collision_one* og_collision = new nu_nu_collision_one(eps, p1);
        og_collision->whole_integral(dens, neutrino, nu_nu_int);
        delete og_collision;
    }
    */
    
    /*
    
    double p1_energy = 3.7;
    double p2_energy = eps->get_value(eps->get_len()-1);
    
    
    for(int j=0; j<eps->get_len(); j++){
        double p3_energy = eps->get_value(j);
        double p4_energy = p1_energy+p2_energy-p3_energy;

        int count = eps->index_below_for_interpolation(p4_energy);
        if(count < 0) {
            std::cout << "Warning: count=" << count << ", p4_energy=" << p4_energy << std::endl;
        }

        int back = 2;
        if (count - back < 0)
            back = count;
        
        double interpolated_p0=0;
        
        //if p4 energy falls below the maximum energy
        if(p4_energy < eps->get_value(eps->get_len()-1)){
            if (count == dens->num_bins()-1){
                back = 3;}
            
            double* eps_values = new double[4];
            double* p_values = new double[4];

            for(int k = 0; k < 4; k++)
            {
                eps_values[k] = eps->get_value(count-back+k);
                p_values[k] = std::log(dens->p0(count-back+k, neutrino));

            }
            interpolated_p0 = std::exp(interpolate_one(p4_energy, 4, eps_values, p_values));
            delete[] eps_values;
            delete[] p_values;


        }
        else{
            count = eps->get_len();
            interpolated_p0 = extrapolate_exponential_one(p4_energy, eps->get_value(count-2), eps->get_value(count-1), dens->p0(count-2, neutrino), dens->p0(count-1, neutrino));
        }
    std::cout << p4_energy << ", " << interpolated_p0 << ", " << std::setprecision(16) << complex<double> (1.,0) -  interpolated_p0 << std::endl;

    }
    */
    
    delete eps;
    delete dens;
    delete[] nu_nu_int;
    return 0;
    
}

double interpolate_one(double x, int N, double* x_vals, double* y_vals)
{
    double res = 0;
    double termj = 1;
    for (int j = 0; j < N; j++)
    {
        termj = 1;
        for(int k = 0; k < N; k++)
            if (j!=k)
                termj *= (x - x_vals[k]) / (x_vals[j] - x_vals[k]);
        termj *= y_vals[j];
        res += termj;
    }
    return res;
}


double extrapolate_exponential_one(double x, double x1, double x2, double y1, double y2){
    //note: this assumes x1<x2, so we expect y1>y2 because this is an exponential decay model
    if(y1==y2){
        return y1;
    }
    
    else{
     //model is Ce^(-ax)
      //  if(y1/y2<1){std::cout << "warning: attempting to take log of something less than 1" << std::endl;}
        if(y1/y2 < 1)
            return extrapolate_linear_one(x, x1, x2, y1, y2);
        if(x1-x2==0){std::cout << "warning: attempting to divide by 0" << x << std::endl;}
        double a = -log(y1/y2) / (x1-x2);
        double C = y1 * exp(a * x1);
        return C * exp(-a * x);
        
    }
    
}

double extrapolate_linear_one(double x, double x1, double x2, double y1, double y2){
    if(x2-x1==0){std::cout << "warning: attempting to divide by 0**" << x << std::endl;}
    double slope = (y2-y1)/(x2-x1);
    double Delta = slope * (x - x2);
    double result = 0;
    
    if (Delta > 0){
        if (y2 < 1){
            return y2 + (1 - y2) * std::tanh(Delta/(1-y2));
        }
        else{
            return 1.0;
        }
    }
    else{
        if (y2 > -1){
            return y2 + (y2 + 1) * std::tanh(Delta/(y2+1));
        }
        else{
            return -1.0;
        }
    }
}


double interpolate_two(double x, double x1, double x2, double y1, double y2){
    //if(x2==x1){std::cout << "Warning: attempting to divide by zero, x1==" << x1 << ", x2==" << x2 << std::endl;}
    //return y1 * abs(x1-x)/abs(x2-x1) + y2 * abs(x2-x)/abs(x2-x1);
    return (y2 - y1)/(x2 - x1) * (x-x1) + y1;
}