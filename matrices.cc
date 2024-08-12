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
double extrapolate(double, double, double, double, double);



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

void matrix::convert_p4_to_interpolated_matrix(density* dens, bool neutrino, double p4_energy, int count){
    dummy_vars* eps = dens->get_E();
    
    double interpolated_p0;
    three_vector* interpolated_p = new three_vector();
    three_vector* p1 = new three_vector();
    three_vector* p2 = new three_vector();
    double temp_result;
    
    //if p4 energy falls below the maximum energy
    if(count<dens->num_bins()){
        dens->p0_p(count-1, neutrino, p1);
        dens->p0_p(count, neutrino, p2);
        
        interpolated_p0 = interpolate(p4_energy, eps->get_value(count-1), eps->get_value(count), dens->p0(count-2, neutrino), dens->p0(count-1, neutrino));
        for(int i=0; i<3; i++){
            temp_result = interpolate(p4_energy, eps->get_value(count-1), eps->get_value(count), p1->get_value(i), p2->get_value(i));
            interpolated_p->set_value(i, temp_result);
        }
    }
    
    else{
        //p0
        interpolated_p0 = extrapolate_exponential(p4_energy, eps->get_value(count-2), eps->get_value(count-1), dens->p0(count-2, neutrino), dens->p0(count-1, neutrino));
        
        dens->p0_p(count-2, neutrino, p1);
        dens->p0_p(count-1, neutrino, p2);
        //px,py
        for(int i=0; i<2; i++){
            temp_result = extrapolate_linear(p4_energy, eps->get_value(count-2), eps->get_value(count-1), p1->get_value(i), p2->get_value(i));
            interpolated_p->set_value(i, temp_result);
        }
        temp_result = extrapolate_exponential(p4_energy, eps->get_value(count-2), eps->get_value(count-1), p1->get_value(2), p2->get_value(2));
        interpolated_p->set_value(2, temp_result);
        
    }
    A0 = complex<double> (0.5 * interpolated_p0, 0);
    A->make_complex(interpolated_p);
    A->multiply_by(complex<double>(0.5,0));
    if(std::isnan(real(A0)) != 0){std::cout <<"constatn multiplier for matrix is nan" << std::endl;}
    
    delete p1;
    delete p2;
    delete interpolated_p;
}

void matrix::convert_p4_to_identity_minus_interpolated_matrix(density* dens, bool neutrino, double p4_energy, int count){
    dummy_vars* eps = dens->get_E();
    
    double interpolated_p0;
    three_vector* interpolated_p = new three_vector();
    three_vector* p1 = new three_vector();
    three_vector* p2 = new three_vector();

    double temp_result;
    
    //if p4 energy falls below the maximum energy: we will do interpolation
    if(count<dens->num_bins()){
        dens->p0_p(count-1, neutrino, p1);
        dens->p0_p(count, neutrino, p2);
        
        interpolated_p0 = interpolate(p4_energy, eps->get_value(count-1), eps->get_value(count), dens->p0(count-2, neutrino), dens->p0(count-1, neutrino));
        for(int i=0; i<3; i++){
            temp_result = interpolate(p4_energy, eps->get_value(count-1), eps->get_value(count), p1->get_value(i), p2->get_value(i));
            interpolated_p->set_value(i, temp_result);
        }
        
    }
    //if p4 is past the end of the array: we will do extrapolation
    else{
        //p0
        interpolated_p0 = extrapolate_exponential(p4_energy, eps->get_value(count-2), eps->get_value(count-1), dens->p0(count-2, neutrino), dens->p0(count-1, neutrino));
        dens->p0_p(count-2, neutrino, p1);
        dens->p0_p(count-1, neutrino, p2);
        
        //px,py
        for(int i=0; i<2; i++){
            temp_result = extrapolate_linear(p4_energy, eps->get_value(count-2), eps->get_value(count-1), p1->get_value(i), p2->get_value(i));
            interpolated_p->set_value(i, temp_result);
        }
        //pz
        temp_result = extrapolate_exponential(p4_energy, eps->get_value(count-2), eps->get_value(count-1), p1->get_value(2), p2->get_value(2));
        interpolated_p->set_value(2, temp_result);
        
    }
    A0 = complex<double> (1 - 0.5 * interpolated_p0, 0);
    if(std::isnan(real(A0)) != 0){std::cout <<"constatn multiplier for matrix is nan" << std::endl;}
    A->make_complex(interpolated_p);
    A->multiply_by(complex<double>(-0.5,0));
    
    delete p1;
    delete p2;
    delete interpolated_p;
}

double interpolate(double x, double x1, double x2, double y1, double y2){
    if(x2==x1){std::cout << "Warning: attempting to divide by zero" << std::endl;}
    return y1 * abs(x1-x)/abs(x2-x1) + y2 * abs(x2-x)/abs(x2-x1);
}

double extrapolate_exponential(double x, double x1, double x2, double y1, double y2){
    //note: this assumes x1<x2, so we expect y1>y2 because this is an exponential decay model
    if(y1==y2){
        return y1;
    }
    
    else{
     //model is Ce^(-ax)
        if(y1/y2<1){std::cout << "warning: attempting to take log of something less than 1" << std::endl;}
        if(x1-x2==0){std::cout << "warning: attempting to divide by 0" << std::endl;}
        double a = -log(y1/y2) / (x1-x2);
        double C = y1 * exp(-a * x1);
        
        return C * exp(-a * x);
        
    }
    
}

double extrapolate_linear(double x, double x1, double x2, double y1, double y2){
    if(x2-x1==0){std::cout << "warning: attempting to divide by 0" << std::endl;}
    double slope = (y2-y1)/(x2-x1);
    return slope * (x - x2) + y2;
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