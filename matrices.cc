#include "matrices.hh"
#include "arrays.hh"
#include "QKE_methods.hh"
#include <complex>
#include <iostream>


using std::cout;
using std::endl;
using std::complex;



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
    linspace_and_gl* eps = dens->get_E();
    
    double interpolated_p0;
    three_vector* interpolated_p = new three_vector();
    three_vector* p1 = new three_vector();
    three_vector* p2 = new three_vector();
    
    //if p4 energy falls below the maximum energy
    if(count<dens->num_bins()){
        double ratio1 = (p4_energy-eps->get_value(count-1))/(eps->get_value(count)-eps->get_value(count-1));
        double ratio2 = (eps->get_value(count)-p4_energy)/(eps->get_value(count)-eps->get_value(count-1));
        
        interpolated_p0 = dens->p0(count-1, neutrino) * ratio1 + dens->p0(count, neutrino) * ratio2;
        
        dens->p0_p(count-1, neutrino, p1);
        dens->p0_p(count, neutrino, p2);
        p1->multiply_by(ratio1);
        p2->multiply_by(ratio2);
        p3->add(p1, p2);
    }
    
    else{
        //assume density is a function of the form Ce^(-ax); we will use the last two points in eps to find C and a
        //given two points (x1,y1) and (x2,y2) on this curve we have a =log(y1/y2)/(x2-x1) and C = y1 * e^(a*x1)
        
        double a = log((dens->p0(count-1, neutrino)-dens->p0(count, neutrino))/(eps->get_value(count)-eps->get_value(count-1)));
        double C = dens->p0(count, neutrino) * exp(a * eps->get_value(count));
        interpolated_p0 = C * exp(-a * p4_energy);
        
        dens->p0_p(count-1, neutrino, p1);
        dens->p0_p(count, neutrino, p2);
        for(int i=0; i<3; i++){
            a = log((p1->get_value(i) - p2->get_value(i))/(eps->get_value(count)-eps->get_value(count-1)));
            C = p2->get_value(i) * exp(a * eps->get_value(count));
            interpolated_p->set_value(i, C * exp(-a * p4_energy));
        }
        
    }
    A0 = complex<double> (0.5 * interpolated_p0, 0);
    A->make-complex(interpolated_p);
    A->multiply_by(complex<double>(0.5,0));
    
    delete p1;
    delete p2;
    delete interpolated_p;
}

void matrix::convert_p4_to_identity_minus_interpolated_matrix(density* dens, bool neutrino, double p4_energy, int count){
    linspace_and_gl* eps = dens->get_E();
    
    double interpolated_p0;
    three_vector* interpolated_p = new three_vector();
    three_vector* p1 = new three_vector();
    three_vector* p2 = new three_vector();
    
    //if p4 energy falls below the maximum energy
    if(count<dens->num_bins()){
        double ratio1 = (p4_energy-eps->get_value(count-1))/(eps->get_value(count)-eps->get_value(count-1));
        double ratio2 = (eps->get_value(count)-p4_energy)/(eps->get_value(count)-eps->get_value(count-1));
        
        interpolated_p0 = dens->p0(count-1, neutrino) * ratio1 + dens->p0(count, neutrino) * ratio2;
        
        dens->p0_p(count-1, neutrino, p1);
        dens->p0_p(count, neutrino, p2);
        p1->multiply_by(ratio1);
        p2->multiply_by(ratio2);
        p3->add(p1, p2);
    }
    
    else{
        //assume density is a function of the form Ce^(-ax); we will use the last two points in eps to find C and a
        //given two points (x1,y1) and (x2,y2) on this curve we have a =log(y1/y2)/(x2-x1) and C = y1 * e^(a*x1)
        
        double a = log((dens->p0(count-1, neutrino)-dens->p0(count, neutrino))/(eps->get_value(count)-eps->get_value(count-1)));
        double C = dens->p0(count, neutrino) * exp(a * eps->get_value(count));
        interpolated_p0 = C * exp(-a * p4_energy);
        
        dens->p_vector(count-1, neutrino, p1);
        dens->p_vector(count, neutrino, p2);
        for(int i=0; i<3; i++){
            a = log((p1->get_value(i) - p2->get_value(i))/(eps->get_value(count)-eps->get_value(count-1)));
            C = p2->get_value(i) * exp(a * eps->get_value(count));
            interpolated_p->set_value(i, C * exp(-a * p4_energy));
        }
        
    }
    A0 = complex<double> (1 - 0.5 * interpolated_p0, 0);
    A->make-complex(interpolated_p);
    A->multiply_by(complex<double>(-0.5,0));
    
    delete p1;
    delete p2;
    delete interpolated_p;
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