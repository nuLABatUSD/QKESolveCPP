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
void all_F_for_p1(density*, bool, int, double***);
double integral_one(density*, bool, int, double**);
double integral_two(density*, bool, int, double**);
double integral_three(density*, bool, int, double**);
double integral_four(density*, bool, int, double**);
double integral_five(density*, bool, int, double**);
double integral_six(density*, bool, int, double**);
double whole_integral(density*, bool, int, double**);

int main(){
    
    
    linspace_for_trap* et = new linspace_for_trap(0.,20, 201);
    double eta_e = 0.2;
    double eta_mu = -0.02;
    density* den = new density(et, eta_e, eta_mu);
    
    
    double*** F_values = new double**[4];
    for(int i=0; i<4; i++){
        F_values[i] = new double*[den->num_bins()];
        for(int j=0; j<den->num_bins(); j++){
            F_values[i][j] = new double[den->num_bins()];   
        }
    }
    
    
    all_F_for_p1(den, false, 5, F_values);
    
    cout << whole_integral(den, false, 5, F_values[0]);
    
    for(int i=0; i<4; i++){
        for(int j=0; j<den->num_bins(); j++){
            delete[] F_values[i][j];
        }
        delete[] F_values[i];
    }
    delete[] F_values;
    
    
    delete et;
    delete den;
    
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

//Note: this function has no safeguard for negative p4, will produce a result with no error
void Fvvsc_components(density* dens, bool neutrino, int p1, int p2, int p3, double* F03, three_vector* F3){
    
    double F01;
    three_vector* F1 = new three_vector();
    double F02;
    three_vector* F2 = new three_vector();
    
    
    Fvvsc_components_term_1(dens, neutrino, p1, p2, p3, &F01, F1);
    Fvvsc_components_term_2(dens, neutrino, p1, p2, p3, &F02, F2);
    
    F2->multiply_by(-1);
    F3->add(F1, F2);
    
    *F03 = F01 - F02;
    
    delete F1;
    delete F2;
}

/*
populates the double*** object F_vals, which is a 4 by eps->N by eps->N matrix
F0: index is [0][p2][p3]
Fx: index is [1][p2][p3]
Fy: index is [2][p2][p3]
Fz: index is [3][p2][p3]
*/

void all_F_for_p1(density* dens, bool neutrino, int p1, double*** F_vals){
    dummy_vars* eps = dens->get_E();
    
    double F0 = 0;
    three_vector* Fxyz = new three_vector();
    
    for(int p2=0; p2<eps->N; p2++){
        for(int p3=0; p3<eps->N; p3++){
            
            if(p1+p2-p3>=0){
                Fvvsc_components(dens, neutrino, p1, p2, p3, &F0, Fxyz);
                F_vals[0][p2][p3] = F0;
                F_vals[1][p2][p3] = Fxyz->get_value(0);
                F_vals[2][p2][p3] = Fxyz->get_value(1);
                F_vals[3][p2][p3] = Fxyz->get_value(2);
            }
        }
    }  
    
    delete Fxyz;
}


double J1(double p1, double p2, double p3){
    return 16./15 * pow(p3, 3) * (10 * pow(p1+p2, 2) - 15 * (p1+p2) * p3 + 6*pow(p3, 2));  
}

double J2(double p1, double p2){
    return 16./15 * pow(p2, 3) * (10 * pow(p1, 2) + 5 * p1 * p2 + pow(p2, 2));  
}

double J3(double p1, double p2, double p3){
    return 16./15 * (pow(p1+p2, 5) - 10 * pow(p1+p2, 2) * pow(p3, 3) + 15 * (p1+p2) * pow(p3, 4) - 6 * pow(p3,5));
}


double integral_one(density* dens, bool neutrino, int p1, double** F_vals){
    dummy_vars* eps = dens->get_E();
    linspace_for_trap* p_2 = new linspace_for_trap(0, eps->get_value(p1), p1+1);
    double p_1_energy = eps->get_value(p1);
    
    dep_vars* dummy_p_2 = new dep_vars(p1+1);
    
    //NOTE: p_2 HAS LENGTH P1+1 EVEN THOUGH WE INTEGRATE OVER P1 POINTS BECAUSE WE ARE REALLY SUPPOSED TO BE INTEGRATING STARTING AT 0 SO IT SHOULD BE P1 POINTS TO MAKE WEIGHTS RIGHT, WE JUST START AT 1 NOT 0 BC INTERIOR INTEGRAL IS 0 WHEN P2=0
    for(int p2=1; p2<=p1; p2++){
        linspace_for_trap* I1_domain = new linspace_for_trap(0, eps->get_value(p2), p2+1);
 
        dep_vars* dummy_p_3 = new dep_vars(p2+1);
        
        for(int p3=0; p3<=p2; p3++){
            dummy_p_3->set_value(p3, F_vals[p2][p3] * J1(p_1_energy, eps->get_value(p2), eps->get_value(p3)));
        }
        dummy_p_2->set_value(p2, I1_domain->integrate(dummy_p_3));
        
        delete dummy_p_3;
        delete I1_domain;
        
    }
    
    double result = p_2->integrate(dummy_p_2);
    
    delete dummy_p_2;
    delete p_2;
    
    return result;
}


//input double** F_vals is 201 by 201 matrix of p2/p3 F values for a particular one of F0, Fx, Fy, Fz
double integral_two(density* dens, bool neutrino, int p1, double** F_vals){
    dummy_vars* eps = dens->get_E();
    double p_1_energy = eps->get_value(p1);
    linspace_for_trap* p_2 = new linspace_for_trap(0, p_1_energy, p1+1);
    double max_energy = eps->get_value(eps->N-1);
    
    dep_vars* dummy_p_2 = new dep_vars(p1+1);
    
    //NOTE: P2 GOES UP TO P1 NOT INCLUDING P1 BC FOR P2=P1 INNER INTEGRAL WILL BE 0
    //ALSO NOTE: EVEN THOUGH THERE ARE P1 POSSIBILITIES FOR P2, WE NEED P1+1 SPACES IN P_2 LINSPACE FOR TRAP ARRAY BC THIS MAKES THE WEIGHTS RIGHT--THE WEIGHTS HAVE TO MATCH THE SPACE WE SHOULD BE INTEGRATING OVER, WHICH INCLUDES P1 AND THEREFORE HAS P1+1 THINGS IN IT
    for(int p2=0; p2<p1; p2++){
        
        linspace_for_trap* I2_domain = new linspace_for_trap(eps->get_value(p2), p_1_energy, p1-p2+1);

        dep_vars* dummy_p_3 = new dep_vars(p1-p2+1);
        
        //NOTE: INTEGRAL IS SET TO 0 FOR ALL P3 GREATER THAN PMAX
        for(int p3=p2; (p3<=p1); p3++){
            dummy_p_3->set_value(p3-p2, F_vals[p2][p3] * J2(p_1_energy, eps->get_value(p2)));
        }
        
        dummy_p_2->set_value(p2, I2_domain->integrate(dummy_p_3));
        
        delete dummy_p_3;
        delete I2_domain;
    }
    
    double result = p_2->integrate(dummy_p_2);
    
    delete dummy_p_2;
    delete p_2;
    
    return result;
}

double integral_three(density* dens, bool neutrino, int p1, double** F_vals){
    dummy_vars* eps = dens->get_E();
    double p_1_energy = eps->get_value(p1);
    linspace_for_trap* p_2 = new linspace_for_trap(0, p_1_energy, p1+1);
    double max_energy = eps->get_value(eps->N-1);
    
    dep_vars* dummy_p_2 = new dep_vars(p1+1);
    
    //NOTE: P2 STARTS AT 1 NOT 0 BC WHEN P2=0, INTERIOR INTEGRAL IS 0
    for(int p2=1; p2<=p1; p2++){
        
        if(eps->get_value(p2)+p_1_energy <= max_energy){
            linspace_for_trap* I3_domain = new linspace_for_trap(p_1_energy, eps->get_value(p2)+p_1_energy, p2+1);
            dep_vars* dummy_p_3 = new dep_vars(p2+1);

            //NOTE: INTEGRAL IS SET TO 0 FOR ALL P3 GREATER THAN PMAX
            for(int p3=p1; (p3<=p1+p2) and (p3<=eps->N-1); p3++){
                dummy_p_3->set_value(p3-p1, F_vals[p2][p3] * J3(p_1_energy, eps->get_value(p2), eps->get_value(p3)));
            }

            dummy_p_2->set_value(p2, I3_domain->integrate(dummy_p_3));

            delete dummy_p_3;
            delete I3_domain;
        }
        
        else{
            
            linspace_for_trap* I3_domain = new linspace_for_trap(p_1_energy, max_energy, eps->N-p2);
            dep_vars* dummy_p_3 = new dep_vars(eps->N-p2);

            //NOTE: INTEGRAL IS SET TO 0 FOR ALL P3 GREATER THAN PMAX
            for(int p3=p1; (p3<=p1+p2) and (p3<=eps->N-1); p3++){
                dummy_p_3->set_value(p3-p1, F_vals[p2][p3] * J3(p_1_energy, eps->get_value(p2), eps->get_value(p3)));
            }

            dummy_p_2->set_value(p2, I3_domain->integrate(dummy_p_3));

            delete dummy_p_3;
            delete I3_domain;
            
        }
    }
    
    double result = p_2->integrate(dummy_p_2);
    
    delete dummy_p_2;
    delete p_2;
    
    return result;
}

double integral_four(density* dens, bool neutrino, int p1, double** F_vals){
    dummy_vars* eps = dens->get_E();
    double p_1_energy = eps->get_value(p1);
    double max_energy = eps->get_value(eps->N-1);
    //NOTE: LENGTH OF INTEGRATING SPACE DEPENDS ON WHAT MAX ENERGY IS, CURRENTLY INTEGRAL SET TO 0 ANYWHERE BEYOND MAX
    linspace_for_trap* p_2 = new linspace_for_trap(p_1_energy, max_energy, eps->N-p1);

    dep_vars* dummy_p_2 = new dep_vars(eps->N-p1);
    
    //NOTE: IF P1=0 THERE WILL BE AN ERROR BECAUSE INNER INTEGRAL WILL BE 0
    for(int p2=p1; p2<=eps->N-1; p2++){
        
        linspace_for_trap* I4_domain = new linspace_for_trap(0, p_1_energy, p1+1);
        dep_vars* dummy_p_3 = new dep_vars(p1+1);
        
        for(int p3=0; p3<=p1; p3++){
            dummy_p_3->set_value(p3, F_vals[p2][p3] * J1(p_1_energy, eps->get_value(p2), eps->get_value(p3)));
        }
        
        dummy_p_2->set_value(p2-p1, I4_domain->integrate(dummy_p_3));

        delete dummy_p_3;
        delete I4_domain;
    }
    
    double result = p_2->integrate(dummy_p_2);
    
    delete dummy_p_2;
    delete p_2;
    
    return result;
}

double integral_five(density* dens, bool neutrino, int p1, double** F_vals){
    dummy_vars* eps = dens->get_E();
    double p_1_energy = eps->get_value(p1);
    double max_energy = eps->get_value(eps->N-1);
    //NOTE: LENGTH OF INTEGRATING SPACE DEPENDS ON WHAT MAX ENERGY IS, CURRENTLY INTEGRAL SET TO 0 ANYWHERE BEYOND MAX
    linspace_for_trap* p_2 = new linspace_for_trap(p_1_energy, max_energy, eps->N-p1);

    dep_vars* dummy_p_2 = new dep_vars(eps->N-p1);
    
    //NOTE: INTEGRAL STARTS AT P1+1 NOT P1 BECAUSE IF P2=P1 THE INTERIOR INTEGRAL IS 0
    for(int p2=p1+1; p2<=eps->N-1; p2++){
        
        linspace_for_trap* I5_domain = new linspace_for_trap(p_1_energy, eps->get_value(p2), p2-p1+1);
        dep_vars* dummy_p_3 = new dep_vars(p2-p1+1);
        
        for(int p3=p1; p3<=p2; p3++){
            dummy_p_3->set_value(p3-p1, F_vals[p2][p3] * J2(eps->get_value(p2), p_1_energy));
           
        }
        
        dummy_p_2->set_value(p2-p1, I5_domain->integrate(dummy_p_3));

        delete dummy_p_3;
        delete I5_domain;
    }
    
    double result = p_2->integrate(dummy_p_2);
    
    delete dummy_p_2;
    delete p_2;
    
    return result;
}

double integral_six(density* dens, bool neutrino, int p1, double** F_vals){
    dummy_vars* eps = dens->get_E();
    double p_1_energy = eps->get_value(p1);
    double max_energy = eps->get_value(eps->N-1);
    //NOTE: LENGTH OF INTEGRATING SPACE DEPENDS ON WHAT MAX ENERGY IS, CURRENTLY INTEGRAL SET TO 0 ANYWHERE BEYOND MAX
    linspace_for_trap* p_2 = new linspace_for_trap(p_1_energy, max_energy, eps->N-p1);

    dep_vars* dummy_p_2 = new dep_vars(eps->N-p1);
    
    //NOTE: IF P1 IS 0 THERE WILL BE AN ERROR BECAUSE THE INNER INTEGRAL IS 0
    //NOTE: INTEGRAL IS NOT INCLUSIVE OF MAX ENERGY BECAUSE IF P2=MAX ENERGY INNTER INTEGRAL IS 0
    for(int p2=p1; p2<eps->N-1; p2++){
        
        //THIS HANDLES SETTING UP INTEGRATION DEPENDING ON IF P1+P1>E_MAX OR NOT
        if(eps->get_value(p2)+p_1_energy <= max_energy){
            linspace_for_trap* I6_domain = new linspace_for_trap(eps->get_value(p2), eps->get_value(p2)+p_1_energy, p1+1);
            dep_vars* dummy_p_3 = new dep_vars(p1+1);
            
            
            for(int p3=p2; (p3<=p2+p1) and (p3<=eps->N-1); p3++){
                    dummy_p_3->set_value(p3-p1, F_vals[p2][p3] * J3(p_1_energy, eps->get_value(p2), eps->get_value(p3)));
                    
                }

                dummy_p_2->set_value(p2-p1, I6_domain->integrate(dummy_p_3));

                delete dummy_p_3;
                delete I6_domain;
        }
        else{
            
            linspace_for_trap* I6_domain = new linspace_for_trap(eps->get_value(p2), max_energy, eps->N-p2);
            dep_vars* dummy_p_3 = new dep_vars(eps->N-p2);
            
            for(int p3=p2; (p3<=p2+p1) and (p3<=eps->N-1); p3++){
                dummy_p_3->set_value(p3-p1, F_vals[p2][p3] * J3(p_1_energy, eps->get_value(p2), eps->get_value(p3)));
                
            }
            
            dummy_p_2->set_value(p2-p1, I6_domain->integrate(dummy_p_3));

            delete dummy_p_3;
            delete I6_domain;
       }
        
        
    }
    
    double result = p_2->integrate(dummy_p_2);
    
    delete dummy_p_2;
    delete p_2;
    
    return result;
}

double whole_integral(density* dens, bool neutrino, int p1, double** F_vals){
    double I = 0;
    dummy_vars* eps = dens->get_E();
    double p_1_energy = eps->get_value(p1);
    
    I += integral_one(dens, neutrino, p1, F_vals);
    I += integral_two(dens, neutrino, p1, F_vals);
    I += integral_three(dens, neutrino, p1, F_vals);
    I += integral_four(dens, neutrino, p1, F_vals);
    I += integral_five(dens, neutrino, p1, F_vals);
    I += integral_six(dens, neutrino, p1, F_vals);
    
    I *= pow(_GF_,2) / (pow(2*_PI_,3) * pow(p_1_energy,2));
    
    
    return I;
    
    
    
}