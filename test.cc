//#include "matrices.hh"
#include <iostream>
#include <complex> 
#include "matrices.hh"
#include "constants.hh"

using std::cout;
using std::endl;
using std::complex;

void Fvvsc_components_term_2(density*, bool, int, int, int, int, double*, three_vector*);
void Fvvsc_components(density*, bool, int, int, int, int, double*, three_vector*);

int main(){
    
    
    linspace_for_trap* et = new linspace_for_trap(0.,20, 201);
    double eta_e = 0.2;
    double eta_mu = -0.02;
    density* den = new density(et, eta_e, eta_mu);
    
    double* F_03;
    three_vector* F3 = new three_vector();
    Fvvsc_components(den, false, 12, 15, 25, 35, F_03, F3);
    
    cout << "Results: " << endl;
    cout << "F0: " << *F_03 << endl << "F: ";
    
    F3->print_all();

    
    delete et;
    delete den;
    delete F3;
    
    return 0;
}

//inputs are density matrix, boolean indicating neutrino/antineutrino, four energy indices, double that will be modified and threevector that will be modified
void Fvvsc_components_term_1(density* dens, bool neutrino, int p1, int p2, int p3, int p4, double* F0, three_vector* F){
    
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();
    matrix* p_4 = new matrix();
    p_1->convert_p_to_identity_minus_matrix(dens, neutrino, p1);
    p_2->convert_p_to_identity_minus_matrix(dens, neutrino, p2);
    p_3->convert_p_to_matrix(dens, neutrino, p3);
    p_4->convert_p_to_matrix(dens, neutrino, p4);
    
    matrix* F_dummy1 = new matrix();
    F_dummy1->matrix_multiply(p_1, p_3);
    matrix* F_dummy2 = new matrix();
    F_dummy2->matrix_multiply(p_2, p_4);

    complex_three_vector* iden = new complex_three_vector();
    matrix* id = new matrix(complex<double> (1,0), iden);
    
    id->multiply_by(F_dummy2->get_A0()*(complex<double> (2,0)));
    
    matrix* F_dummy3 = new matrix();
    F_dummy3->matrix_add(F_dummy2, id);
    
    matrix* F_dummy4 = new matrix();
    F_dummy4->matrix_multiply(F_dummy1, F_dummy3);
    
    complex<double> comp_F0 = F_dummy4->get_A0();
    complex_three_vector* comp_F = F_dummy4->get_A();
    
    comp_F->multiply_by(2);
    
    *F0 = 2*real(comp_F0);
    F->make_real(comp_F);
    
    
    delete F_dummy1;
    delete F_dummy2;
    delete id;
    delete iden;
    delete F_dummy3;
    delete F_dummy4;
    delete p_1;
    delete p_2;
    delete p_3;
    delete p_4;
}


void Fvvsc_components_term_2(density* dens, bool neutrino, int p1, int p2, int p3, int p4, double* F0, three_vector* F){
    
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();
    matrix* p_4 = new matrix();
    p_1->convert_p_to_matrix(dens, neutrino, p1);
    p_2->convert_p_to_matrix(dens, neutrino, p2);
    p_3->convert_p_to_identity_minus_matrix(dens, neutrino, p3);
    p_4->convert_p_to_identity_minus_matrix(dens, neutrino, p4);
    
    matrix* F_dummy1 = new matrix();
    F_dummy1->matrix_multiply(p_1, p_3);
    matrix* F_dummy2 = new matrix();
    F_dummy2->matrix_multiply(p_2, p_4);

    complex_three_vector* iden = new complex_three_vector();
    matrix* id = new matrix(complex<double> (1,0), iden);
    
    id->multiply_by(F_dummy2->get_A0()*(complex<double> (2,0)));
    
    matrix* F_dummy3 = new matrix();
    F_dummy3->matrix_add(F_dummy2, id);
    
    matrix* F_dummy4 = new matrix();
    F_dummy4->matrix_multiply(F_dummy1, F_dummy3);
    
    complex<double> comp_F0 = F_dummy4->get_A0();
    complex_three_vector* comp_F = F_dummy4->get_A();
    
    comp_F->multiply_by(2);
    
    *F0 = 2*real(comp_F0);
    F->make_real(comp_F);
    
    
    delete F_dummy1;
    delete F_dummy2;
    delete id;
    delete iden;
    delete F_dummy3;
    delete F_dummy4;
    delete p_1;
    delete p_2;
    delete p_3;
    delete p_4;
}

void Fvvsc_components(density* dens, bool neutrino, int p1, int p2, int p3, int p4, double* F03, three_vector* F3){
    
    double F01 = 0;
    three_vector* F1 = new three_vector();
    double F02 = 0;
    three_vector* F2 = new three_vector();
    
    Fvvsc_components_term_2(dens, neutrino, p1, p2, p3, p4, &F02, F2);
    Fvvsc_components_term_1(dens, neutrino, p1, p2, p3, p4, &F01, F1);
    
    F2->multiply_by(-1);
    F3->add(F1, F2);
    *F03 = F01 - F02;
    
    delete F1;
    delete F2;
}


void integral_one(density* dens, bool neutrino, int p1, double* sumF0, double* sumFx, double* sumFy, double* sumFz){
    dummy_vars* eps = dens->get_E();
    double p1_energy = eps->get_value(p1);
    
    double* F01;
    
    three_vector* F1 = new three_vector();
    
    for(int p2=0; p2<p1; p2++){
        
        for(int p3=0; p3<p2; p3++){
            Fvvsc_components(dens, neutrino, p1, p2, p3, p1+p2-p3, F01, F1);
            
            double p3_energy = eps->get_value(p3);
            *sumF0 += *F01 + 16./15. * pow(p3_energy,3) * (10 * pow(p1_energy,2) - 5 * p1_energy * p3_energy + pow(p3_energy,2));
            *sumFx += F1->get_value(0) + 16./15. * pow(p3_energy,3) * (10 * pow(p1_energy,2) - 5 * p1_energy * p3_energy + pow(p3_energy,2));
            *sumFy += F1->get_value(1) + 16./15. * pow(p3_energy,3) * (10 * pow(p1_energy,2) - 5 * p1_energy * p3_energy + pow(p3_energy,2));
            *sumFz += F1->get_value(2) + 16./15. * pow(p3_energy,3) * (10 * pow(p1_energy,2) - 5 * p1_energy * p3_energy + pow(p3_energy,2));
        }
        
        
        for(int p3=p2; p3<p1; p3++){
            Fvvsc_components(dens, neutrino, p1, p2, p3, p1+p2-p3, F01, F1);
            
            double p2_energy = eps->get_value(p2);
            double p3_energy = eps->get_value(p3);
            
            *sumF0 += *F01 + 16./15. * pow(p2_energy,2) * (10 * pow(p1_energy - p2_energy,2) + 15 * (p1_energy - p3_energy) * p2_energy + 6 * pow(p2_energy,2));
            *sumFx += F1->get_value(0) + 16./15. * pow(p2_energy,2) * (10 * pow(p1_energy - p2_energy,2) + 15 * (p1_energy - p3_energy) * p2_energy + 6 * pow(p2_energy,2));
            *sumFy += F1->get_value(1) + 16./15. * pow(p2_energy,2) * (10 * pow(p1_energy - p2_energy,2) + 15 * (p1_energy - p3_energy) * p2_energy + 6 * pow(p2_energy,2));
            *sumFz += F1->get_value(2) + 16./15. * pow(p2_energy,2) * (10 * pow(p1_energy - p2_energy,2) + 15 * (p1_energy - p3_energy) * p2_energy + 6 * pow(p2_energy,2));
            
        }
        
        for(int p3=p1; p3<p1+p2; p3++){
            Fvvsc_components(dens, neutrino, p1, p2, p3, p1+p2-p3, F01, F1);
            
            double p2_energy = eps->get_value(p2);
            double p3_energy = eps->get_value(p3);
            
            *sumF0 += *F01 + 16./15. * (pow(p1_energy - p3_energy,5) + 10 * pow(p1_energy - p3_energy,2) * pow(p2_energy,3) + 15 * (p1_energy - p3_energy) * pow(p2_energy,4) + 6 * pow(p2_energy,5));
            *sumFx += F1->get_value(0) + 16./15. * (pow(p1_energy - p3_energy,5) + 10 * pow(p1_energy - p3_energy,2) * pow(p2_energy,3) + 15 * (p1_energy - p3_energy) * pow(p2_energy,4) + 6 * pow(p2_energy,5));
            *sumFy += F1->get_value(1) + 16./15. * (pow(p1_energy - p3_energy,5) + 10 * pow(p1_energy - p3_energy,2) * pow(p2_energy,3) + 15 * (p1_energy - p3_energy) * pow(p2_energy,4) + 6 * pow(p2_energy,5));
            *sumFz += F1->get_value(2) + 16./15. * (pow(p1_energy - p3_energy,5) + 10 * pow(p1_energy - p3_energy,2) * pow(p2_energy,3) + 15 * (p1_energy - p3_energy) * pow(p2_energy,4) + 6 * pow(p2_energy,5));
        }
    }
    delete F1;
    
    *sumF0 *= pow(_GF_,2) / (4 * pow(2*_PI_, 3) * pow(p1_energy,2));
    *sumFx *= pow(_GF_,2) / (4 * pow(2*_PI_, 3) * pow(p1_energy,2));
    *sumFy *= pow(_GF_,2) / (4 * pow(2*_PI_, 3) * pow(p1_energy,2));
    *sumFz *= pow(_GF_,2) / (4 * pow(2*_PI_, 3) * pow(p1_energy,2));

}