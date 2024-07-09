//#include "matrices.hh"
#include <iostream>
#include <complex> 
#include "matrices.hh"
#include "constants.hh"
#include <iomanip>
#include "gl_vals.hh"

using std::cout;
using std::endl;
using std::complex;

void Fvvsc_components_term_1(density*, bool, int, int, int, double*, three_vector*);
void Fvvsc_components_term_2(density*, bool, int, int, int, double*, three_vector*);
void Fvvsc_components(density*, bool, int, int, int, double*, three_vector*);
void all_F_for_p1(density*, bool, int, double***);
double interior_integral(density*, bool, int, int, double**);
double exterior_integral(density*, bool, int, double**);

int main(){
    
    
    linspace_for_trap* et = new linspace_for_trap(0.,20, 201);
    double eta_e = 0.2;
    double eta_mu = -0.02;
    linspace_and_gl* new_et = new linspace_and_gl(0,20,201,5);
    
    
    density* new_den = new density(new_et, eta_e, eta_mu);
    density* den = new density(et, eta_e, eta_mu);
    
    double*** F_values = new double**[4];
    for(int i=0; i<4; i++){
        F_values[i] = new double*[new_den->num_bins()];
        for(int j=0; j<den->num_bins(); j++){
            F_values[i][j] = new double[den->num_bins()];   
        }
    }
    
    int p1 = 200;
    all_F_for_p1(den, true, p1, F_values);
    
    cout << "og version: " << exterior_integral(den, true, p1, F_values[0]) << endl;
    
    for(int i=0; i<4; i++){
        for(int j=0; j<den->num_bins(); j++){
            delete[] F_values[i][j];
        }
        delete[] F_values[i];
    }
    delete[] F_values;
    
    
    delete et;
    delete new_et;
    delete new_den;
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


double interior_integral(density* dens, bool neutrino, int p1, int p2, double** F_vals){
    
    dummy_vars* eps = dens->get_E();
    double p_1_energy = eps->get_value(p1);
    double max_energy = eps->get_value(eps->N-1);
    
    if(eps->get_value(p2)+p_1_energy <= max_energy){
        
        linspace_for_trap* I_domain = new linspace_for_trap(0, eps->get_value(p2)+p_1_energy, p2+p1+1);
        dep_vars* dummy_p_3 = new dep_vars(p2+p1+1);

        for(int p3=0; p3<p1; p3++){
             dummy_p_3->set_value(p3, F_vals[p2][p3] * J1(p_1_energy, eps->get_value(p2), eps->get_value(p3)));
        }
        
        for(int p3=p1; p3<p2; p3++){
            if(p2<p1){
                dummy_p_3->set_value(p3, F_vals[p2][p3] * J2(p_1_energy, eps->get_value(p2)));
            }
            else{
                dummy_p_3->set_value(p3, F_vals[p2][p3] * J2(eps->get_value(p2), p_1_energy));
            }
        }
        
        for(int p3=p2; p3<=p1+p2; p3++){
            dummy_p_3->set_value(p3, F_vals[p2][p3] * J3(p_1_energy, eps->get_value(p2), eps->get_value(p3)));
        }

        double result = I_domain->integrate(dummy_p_3);

        delete dummy_p_3;
        delete I_domain;
        
        return result;
        
    }
    else{
        
        int count = 0;
        for(int i=0; i<eps->N; i++){
            if(eps->get_value(i)<=eps->get_value(p2)+p_1_energy){
                count++;
            }
        }
        count -= 1;
        cout << "count: " << count;
        //count now gives the index of the last energy value we will take from eps to use as part of I_domain
        
        //I_domain will have count+2 elements because count+1 accounts for all of the elements that will be taken from eps and then one more for p1+p2 itself
        dummy_vars* I_domain = new dummy_vars(count+2);
        for(int i=0; i<=count; i++){
            I_domain->set_value(i, eps->get_value(i));
        }
        I_domain->set_value(count+1,eps->get_value(p2)+p_1_energy);
        I_domain->set_trap_weights();
        
        dep_vars* dummy_p_3 = new dep_vars(count+2);

        for(int p3=0; p3<p1; p3++){
             dummy_p_3->set_value(p3, F_vals[p2][p3] * J1(p_1_energy, eps->get_value(p2), eps->get_value(p3)));
        }

        for(int p3=p1; p3<p2; p3++){
            if(p2<p1){
                dummy_p_3->set_value(p3, F_vals[p2][p3] * J2(p_1_energy, eps->get_value(p2)));
            }
            else{
                dummy_p_3->set_value(p3, F_vals[p2][p3] * J2(eps->get_value(p2), p_1_energy));
            }
        }

        for(int p3=p2; p3<=count; p3++){
            dummy_p_3->set_value(p3, F_vals[p2][p3] * J3(p_1_energy, eps->get_value(p2), eps->get_value(p3)));
        }
        
        //note will throw an error if p1+p2 is past end of eps
        double interpolated_F_val = F_vals[p2][count] * (eps->get_value(p2)+p_1_energy-eps->get_value(count))/(eps->get_value(count+1)-eps->get_value(count)) + F_vals[p2][count+1] * (eps->get_value(count+1)-eps->get_value(p2)-p_1_energy)/(eps->get_value(count+1)-eps->get_value(count));            
      
        
        dummy_p_3->set_value(count+1, interpolated_F_val * J3(p_1_energy, eps->get_value(p2), eps->get_value(p2)+p_1_energy));
        
        double result = eps->integrate(dummy_p_3);

        delete dummy_p_3;
        
        return result;

   }
    
}

double exterior_integral(density* dens, bool neutrino, int p1, double** F_vals){
    if (p1==0){return 0;}
    
    dummy_vars* eps = dens->get_E();
    double p_1_energy = eps->get_value(p1);
    
    dep_vars* dummy_p_2 = new dep_vars(eps->N);

    for(int p2=0; p2<eps->N; p2++){
        dummy_p_2->set_value(p2, interior_integral(dens, neutrino, p1, p2, F_vals));
    }

    double result = eps->integrate(dummy_p_2);
    result *= pow(_GF_,2) / (pow(2*_PI_,3) * pow(p_1_energy,2));
    delete dummy_p_2;
    return result;
}
