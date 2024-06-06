//#include "matrices.hh"
#include <iostream>
#include <complex> 
#include "matrices.hh"
#include "constants.hh"
#include <iomanip>

using std::cout;
using std::endl;
using std::complex;

void Fvvsc_components_term_1(density*, bool, int, int, int, double*, three_vector*);
void Fvvsc_components_term_2(density*, bool, int, int, int, double*, three_vector*);
void Fvvsc_components(density*, bool, int, int, int, double*, three_vector*);

int main(){
    
    
    linspace_for_trap* et = new linspace_for_trap(0.,20, 201);
    double eta_e = 0.2;
    double eta_mu = -0.02;
    density* den = new density(et, eta_e, eta_mu);
    
    double* F_03;
    three_vector* F3 = new three_vector();
    Fvvsc_components(den, false, 1, 2, 2, F_03, F3);
    
    cout << "full thing: " << endl;
    cout << "F0: " << *F_03 << endl << "F: ";    
    F3->print_all();
    
    delete et;
    delete den;
    delete F3;
    
    return 0;
}

//inputs are density matrix, boolean indicating neutrino/antineutrino, four energy indices, double that will be modified and threevector that will be modified
void Fvvsc_components_term_1(density* dens, bool neutrino, int p1, int p2, int p3, double* F0, three_vector* F){
    
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();
    matrix* p_4 = new matrix();
    p_1->convert_p_to_identity_minus_matrix(dens, neutrino, p1);
    p_2->convert_p_to_identity_minus_matrix(dens, neutrino, p2);
    p_3->convert_p_to_matrix(dens, neutrino, p3);
    
    
    //this clause makes the p4 matrix 0 if p4 is bigger than max energy
    if (p1+p2-p3 <= dens->num_bins()){
        p_4->convert_p_to_matrix(dens, neutrino, p1+p2-p3);
    }
    
    /*
    
    p_1 = 1-rho_1
    p_2 = 1-rho_2
    p_3 = rho_3
    p_4 = rho_4
    
    F_dummy1 = (1-rho_2)(rho_4)
    id = tr((1-rho_2)(rho_4)) = 2*A_0*1
    F_dummy2 = F_dummy1 + id
    F_dummy3 = (rho_3) * F_dummy2
    F_dummy4 = (1-rho_1) * F_dummy3
    
    */
    
    
    matrix* F_dummy1 = new matrix();
    F_dummy1->matrix_multiply(p_2, p_4);
    
    matrix* id = new matrix(true);
    id->multiply_by(F_dummy1->get_A0()*(complex<double> (2,0)));
    
    matrix* F_dummy2 = new matrix();
    F_dummy2->matrix_add(F_dummy1, id);
    
    matrix* F_dummy3 = new matrix();
    F_dummy3->matrix_multiply(p_3, F_dummy2);
    
    
    matrix* F_dummy4 = new matrix();
    F_dummy4->matrix_multiply(p_1, F_dummy3);
    
    complex<double> comp_F0 = F_dummy4->get_A0();
    complex_three_vector* comp_F = F_dummy4->get_A();
    
    comp_F->multiply_by(2);
    
    *F0 = 2*real(comp_F0);
    F->make_real(comp_F);
    
    
    delete F_dummy1;
    delete F_dummy2;
    delete id;
    delete F_dummy3;
    delete F_dummy4;
    delete p_1;
    delete p_2;
    delete p_3;
    delete p_4;
    
}


void Fvvsc_components_term_2(density* dens, bool neutrino, int p1, int p2, int p3, double* F0, three_vector* F){
    
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();
    
    matrix* p_4 = new matrix();
    p_1->convert_p_to_matrix(dens, neutrino, p1);
    p_2->convert_p_to_matrix(dens, neutrino, p2);
    p_3->convert_p_to_identity_minus_matrix(dens, neutrino, p3);
    
    //this clause makes the p4 matrix 0 if p4 is bigger than max energy
    
    if (p1+p2-p3 <= dens->num_bins()){
        p_4->convert_p_to_identity_minus_matrix(dens, neutrino, p1+p2-p3);
    }
    
    else{
        matrix* p_4 = new matrix(true);
    }
    
    /*
    p_1 = rho_1
    p_2 = rho_2
    p_3 = 1-rho_3
    p_4 = 1-rho_4
    
    F_dummy1 = (rho_2)(1-rho_4)
    id = tr((rho_2)(1-rho_4)) = 2*A_0*1
    F_dummy2 = F_dummy1 + id
    F_dummy3 = (1-rho_3) * F_dummy2
    F_dummy4 = (rho_1) * F_dummy3
    
    */
    
    matrix* F_dummy1 = new matrix();
    F_dummy1->matrix_multiply(p_2, p_4);

    matrix* id = new matrix(true);
    id->multiply_by(F_dummy1->get_A0()*(complex<double> (2,0)));
    
    matrix* F_dummy2 = new matrix();
    F_dummy2->matrix_add(F_dummy1, id);
    
    matrix* F_dummy3 = new matrix();
    F_dummy3->matrix_multiply(p_3,F_dummy2);

    matrix* F_dummy4 = new matrix();
    F_dummy4->matrix_multiply(p_1,F_dummy3); 

    
    complex<double> comp_F0 = F_dummy4->get_A0();
    complex_three_vector* comp_F = F_dummy4->get_A();
    
    comp_F->multiply_by(2);
    
    *F0 = 2*real(comp_F0);
    F->make_real(comp_F);
    
    
    delete F_dummy1;
    delete F_dummy2;
    delete id;
    delete F_dummy3;
    delete F_dummy4;
    delete p_1;
    delete p_2;
    delete p_3;
    delete p_4;
}

void Fvvsc_components(density* dens, bool neutrino, int p1, int p2, int p3, double* F03, three_vector* F3){
    
    double F01;
    three_vector* F1 = new three_vector();
    double F02;
    three_vector* F2 = new three_vector();
    
    
    Fvvsc_components_term_1(dens, neutrino, p1, p2, p3, &F01, F1);
    Fvvsc_components_term_2(dens, neutrino, p1, p2, p3, &F02, F2);
    
    cout << "term 1: " << endl;
    cout << "F0: " << std::setprecision (60) << F01 << endl << "F: ";    
    F1->print_all();
    cout << "term 2: " << endl;
    cout << "F0: " << F02 << endl << "F: ";    
    F2->print_all();
    
    F2->multiply_by(-1);
    F3->add(F1, F2);
    
    *F03 = F01 - F02;
    
    delete F1;
    delete F2;
}


/* 

F_vals is N*N length dep vars object (this comes from there being N options for p2 and N options for p3 so N^2 total combinations
there is one F_vals object for every p1 value

the dep_vars object is ordered first by p2, then by p3--so it looks like {(val for p1=0,p2=0), (val for p1=0,p2=1), ... (val for p1=0,p2=N), (val for p1=1,p2=0), ... (val for p1=N,p2=N)}
each val in and of itself is four values

so to get the F0 value for p2=i and p3=j find F_vals->get_value(i*N+j)
to get the Fx value find F_vals->get_value(i*N+j+1)
to get the Fy value find F_vals->get_value(i*N+j+2)
to get the Fz value find F_vals->get_value(i*N+j+3)

Note: it is already build into Fvvsc_components to set rho_4 to 0 when p4 is outside max energy

*/

/*
void all_F_for_p1(density* dens, bool neutrino, int p1, dep_vars* F_vals){
    dummy_vars* eps = dens->get_E();
    
    double* F0;
    three_vector* Fxyz = new three_vector();
    
    for(int p2=0; p2<eps->N; p2++){
        for (int p3=0; p3<eps->N; p3++){
            Fvvsc_components(dens, neutrino, p1, p2, p3, F0, Fxyz);

            F_vals->set_value(p2*eps->N + p3, *F0);
            F_vals->set_value(p2*eps->N + p3 + 1, Fxyz->get_value(0));
            F_vals->set_value(p2*eps->N + p3 + 2, Fxyz->get_value(1));
            F_vals->set_value(p2*eps->N + p3 + 3, Fxyz->get_value(2));
        }
        
    }
    
    delete Fxyz;
}
*/

double K1(double p1, double p3){
    return 16./15 * pow(p3, 3) * (10 * pow(p1, 2) - 5 * p1 * p3 + pow(p3, 2));  
}

double K2(double p1, double p2, double p3){
    return 16./15 * pow(p2, 3) * (10 * pow(p1-p3, 2) + 15 * (p1-p3) * p3 + 6*pow(p2, 2));  
}

double K3(double p1, double p2, double p3){
    return 16./15 * (pow(p1-p3, 5) + 10 * pow(p1-p3, 2) * pow(p2, 3) + 15 * (p1-p3) * pow(p2, 4) + 6 * pow(p2,5));
}

/*
double F0_integral(density* dens, bool neutrino, int p1){
    dummy_vars* eps = dens->get_E();
    linspace_for_trap* p_2 = new linspace_for_trap(0, eps->get_value(p1), eps->N);
    
    double* F0;
    three_vector* F = new three_vector();
    double interior_integrals = 0;
    
    linspace_for_trap* I1_domain = new linspace_for_trap(0, eps->get_value(p2), p2+1);
    
    
    dep_vars* dummy_p_3 = new dep_vars(p2+1);
    for(int i=0; i<eps->N; i++){
        dummy_p_3->set_value(i) = F0 * K1(p1,eps->get_value(i));
    }
    interior_integrals += I1_domain->integrate(dummy_p_3);
    
    for(int i=p2; i<p1; i++){
        
    }
    
                                          
    
    
}
*/

/*

void integral_one(density* dens, bool neutrino, int p1, double* sumF0, double* sumFx, double* sumFy, double* sumFz){
    dummy_vars* eps = dens->get_E();
    double p1_energy = eps->get_value(p1);
    
    double* F01;
    
    three_vector* F1 = new three_vector();
    
    for(int p2=0; p2<p1; p2++){
        
        for(int p3=0; p3<p2; p3++){
            //Fvvsc_components(dens, neutrino, p1, p2, p3, F01, F1);
            
            double p3_energy = eps->get_value(p3);
            *sumF0 += *F01 + 16./15. * pow(p3_energy,3) * (10 * pow(p1_energy,2) - 5 * p1_energy * p3_energy + pow(p3_energy,2));
            *sumFx += F1->get_value(0) + 16./15. * pow(p3_energy,3) * (10 * pow(p1_energy,2) - 5 * p1_energy * p3_energy + pow(p3_energy,2));
            *sumFy += F1->get_value(1) + 16./15. * pow(p3_energy,3) * (10 * pow(p1_energy,2) - 5 * p1_energy * p3_energy + pow(p3_energy,2));
            *sumFz += F1->get_value(2) + 16./15. * pow(p3_energy,3) * (10 * pow(p1_energy,2) - 5 * p1_energy * p3_energy + pow(p3_energy,2));
        }
        
        
        for(int p3=p2; p3<p1; p3++){
            //Fvvsc_components(dens, neutrino, p1, p2, p3, F01, F1);
            
            double p2_energy = eps->get_value(p2);
            double p3_energy = eps->get_value(p3);
            
            *sumF0 += *F01 + 16./15. * pow(p2_energy,2) * (10 * pow(p1_energy - p2_energy,2) + 15 * (p1_energy - p3_energy) * p2_energy + 6 * pow(p2_energy,2));
            *sumFx += F1->get_value(0) + 16./15. * pow(p2_energy,2) * (10 * pow(p1_energy - p2_energy,2) + 15 * (p1_energy - p3_energy) * p2_energy + 6 * pow(p2_energy,2));
            *sumFy += F1->get_value(1) + 16./15. * pow(p2_energy,2) * (10 * pow(p1_energy - p2_energy,2) + 15 * (p1_energy - p3_energy) * p2_energy + 6 * pow(p2_energy,2));
            *sumFz += F1->get_value(2) + 16./15. * pow(p2_energy,2) * (10 * pow(p1_energy - p2_energy,2) + 15 * (p1_energy - p3_energy) * p2_energy + 6 * pow(p2_energy,2));
            
        }
        
        for(int p3=p1; p3<p1+p2; p3++){
            //Fvvsc_components(dens, neutrino, p1, p2, p3, F01, F1);
            
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
*/
