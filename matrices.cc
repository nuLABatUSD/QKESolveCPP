#include "matrices.hh"
#include "arrays.hh"
#include "QKE_methods.hh"
#include <complex>
#include <iostream>


using std::cout;
using std::endl;
using std::complex;
using std::abs;

double interpolate(double, double, double, double, double);
double extrapolate_exponential(double, double, double, double, double);
double extrapolate_linear(double, double, double, double, double);

double interpolate(double, int, double*, double*);


matrix::matrix(){
    A0 = complex<double> (0,0);
    A = new complex_three_vector();
    
}

matrix::matrix(complex<double> c, complex_three_vector* C){
    A0 = c;
    A = new complex_three_vector();
    for(int i=0; i<3; i++){
       A->set_value(i,C->get_value(i));
    }
    
}

matrix::matrix(bool id){
    if (id==true){
        A0 = complex<double> (1,0);
        A = new complex_three_vector();
    }
    
}

complex<double> matrix::get_A0(){
    return A0;
}

complex_three_vector* matrix::get_A(){
    return A;
}

void matrix::convert_p_to_matrix(double p0, three_vector* p0p){
    A0 = complex<double> (0.5 * p0,0);
    A->make_complex(p0p);
    A->multiply_by(complex<double>(0.5,0));
}

void matrix::convert_p_to_matrix(density* dens, bool neutrino, int i){
    A0 = complex<double> (0.5 * dens->p0(i, neutrino),0);
    three_vector* p0p = new three_vector();
    dens->p0_p(i, neutrino, p0p);
    A->make_complex(p0p);
    A->multiply_by(complex<double>(0.5,0));
    delete p0p;
}

void matrix::convert_p_to_identity_minus_matrix(double p0, three_vector* p0p){
    A0 = complex<double> (1 - 0.5 * p0,0);
    A->make_complex(p0p);
    A->multiply_by(complex<double>(-0.5,0));
}

void matrix::convert_p_to_identity_minus_matrix(density* dens, bool neutrino, int i){
    A0 = complex<double> (1 - 0.5 * dens->p0(i, neutrino),0);
    three_vector* p0p = new three_vector();
    dens->p0_p(i, neutrino, p0p);
    A->make_complex(p0p);
    A->multiply_by(complex<double>(-0.5,0));
    delete p0p;
}

void matrix::convert_p4_to_interpolated_matrix(density* dens, bool neutrino, double p4_energy){
    dummy_vars* eps = dens->get_E();
    
    double interpolated_p0;
    three_vector* interpolated_p = new three_vector();
    three_vector** p_interp = new three_vector*[4];
    for(int j = 0; j < 4; j++)
        p_interp[j] = new three_vector();
    
    double temp_result;
    
    double eps_values[4];
    double p_values[4];
    
    int count = eps->index_below_for_interpolation(p4_energy);
    if(count < 0) {
        std::cout << "Warning: count=" << count << ", p4_energy=" << p4_energy << std::endl;
    }
    
    int back = 2;
    if (count - back < 0)
        back = count;
    //if p4 energy falls below the maximum energy
    if(p4_energy < eps->get_value(eps->get_len()-1)){
        if (count == dens->num_bins()-1)
            back = 3;
            
        for(int j = 0; j < 4; j++)
        {
            eps_values[j] = eps->get_value(count-back+j);
            p_values[j] = std::log(dens->p0(count-back+j, neutrino));
            
            dens->p0_p(count-back+j, neutrino, p_interp[j]);
        }
        interpolated_p0 = std::exp(interpolate(p4_energy, 4, eps_values, p_values));
        //std::cout << interpolated_p0 << std::endl;
        
        for(int i=0; i<3; i++){
            for(int j = 0; j < 4; j++)
                p_values[j] = p_interp[j]->get_value(i);
                temp_result = interpolate(p4_energy, 4, eps_values, p_values);
                //std::cout << temp_result << ", ";
                interpolated_p->set_value(i, temp_result);
        }
        //std::cout << " " << std::endl;
    }
    
    else{
        //p0
        count = eps->get_len();

        interpolated_p0 = extrapolate_exponential(p4_energy, eps->get_value(count-2), eps->get_value(count-1), dens->p0(count-2, neutrino), dens->p0(count-1, neutrino));
        //std::cout << interpolated_p0 << std::endl;
        
        
        dens->p_vector(count-2, neutrino, p_interp[0]);
        dens->p_vector(count-1, neutrino, p_interp[1]);
        
        //px,py,pz
        for(int i=0; i<3; i++){
            temp_result = extrapolate_linear(p4_energy, eps->get_value(count-2), eps->get_value(count-1), p_interp[0]->get_value(i), p_interp[1]->get_value(i));
            //std::cout << temp_result << ", ";
            
            temp_result *= interpolated_p0;
            interpolated_p->set_value(i, temp_result);
        }
       // std::cout << " " << std::endl;
        
    }
    A0 = complex<double> (0.5 * interpolated_p0, 0);
    A->make_complex(interpolated_p);
    A->multiply_by(complex<double>(0.5,0));
    if(std::isnan(real(A0)) != 0){std::cout <<"constatn multiplier for matrix is nan, count=" << count << std::endl << p4_energy << ", " << count << std::endl;}
    
    
    for(int j = 0; j < 4; j++)
        delete p_interp[j];
    delete[] p_interp;
    delete interpolated_p;
}

void matrix::convert_p4_to_identity_minus_interpolated_matrix(density* dens, bool neutrino, double p4_energy)
{
    convert_p4_to_interpolated_matrix(dens, neutrino, p4_energy);
    
    convert_this_to_identity_minus_this();

}

void matrix::convert_this_to_identity_minus_this()
{
    A0 = complex<double> (1.0) - A0;
    A->multiply_by(complex<double>(-1.,0.));
}


double interpolate(double x, int N, double* x_vals, double* y_vals)
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

double interpolate(double x, double x1, double x2, double y1, double y2){
    //if(x2==x1){std::cout << "Warning: attempting to divide by zero, x1==" << x1 << ", x2==" << x2 << std::endl;}
    return y1 * abs(x1-x)/abs(x2-x1) + y2 * abs(x2-x)/abs(x2-x1);
}

double extrapolate_exponential(double x, double x1, double x2, double y1, double y2){
    //note: this assumes x1<x2, so we expect y1>y2 because this is an exponential decay model
    if(y1==y2){
        return y1;
    }
    
    else{
     //model is Ce^(-ax)
      //  if(y1/y2<1){std::cout << "warning: attempting to take log of something less than 1" << std::endl;}
        if(y1/y2 < 1)
            return extrapolate_linear(x, x1, x2, y1, y2);
        if(x1-x2==0){std::cout << "warning: attempting to divide by 0" << x << std::endl;}
        double a = -log(y1/y2) / (x1-x2);
        double C = y1 * exp(a * x1);
        return C * exp(-a * x);
        
    }
    
}

double extrapolate_linear(double x, double x1, double x2, double y1, double y2){
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

void matrix::print_all(){
    cout << "A_0: " << A0 << endl;
    cout << "A: (" << A->get_value(0) << ", " << A->get_value(1) << ", " << A->get_value(2) << ")" << endl;
}

//add two matrices, modifies whatever matrix it is called on
void matrix::matrix_add(matrix* C1, matrix* C2){
    complex_three_vector* C1_A = C1->get_A();
    complex_three_vector* C2_A = C2->get_A();
    A->add(C1_A, C2_A);
    A0 = C1->get_A0()+C2->get_A0();
    
    
}

//multiply whole matrix by complex number
void matrix::multiply_by(complex<double> c){
    A0 *= c;
    A->multiply_by(c);
}

//multiply two matrices, modifies whatever matrix it is called on
void matrix::matrix_multiply(matrix* C1, matrix* C2){
    complex_three_vector* C1_A = new complex_three_vector(C1->get_A());
    complex_three_vector* C2_A = new complex_three_vector(C2->get_A());
    complex<double> C1_A0 = C1->get_A0();
    complex<double> C2_A0 = C2->get_A0();
    
    complex_three_vector* cp = new complex_three_vector();
    cp->set_cross_product(C1_A, C2_A);
    
    
    A0 = C1_A0*C2_A0 + C1_A->dot_with(C2_A);
    C1_A->multiply_by(C2_A0);
    C2_A->multiply_by(C1_A0);
    A->add(C1_A, C2_A);
    A->add(A,cp);
    
    
    delete cp;
    delete C1_A;
    delete C2_A;
    
}

matrix::~matrix(){
   delete A; 
}