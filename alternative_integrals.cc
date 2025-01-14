#include "alternative_integrals.hh"
#include "QKE_methods.hh"
#include <iostream>
#include "constants.hh"
#include <cmath>
#include "thermodynamics.hh"
#include "gl_vals.hh"
#include "matrices.hh"
#include <fstream>

nu_nu_collision_one::nu_nu_collision_one(linspace_and_gl* e, int p1_index){
    eps = new linspace_and_gl(e);
    p1 = p1_index;
    count = 0;
    outer_vals = new dep_vars(eps->get_len());
    inner_vals = new dep_vars*[eps->get_len()];
    p3_vals = new dummy_vars*[eps->get_len()];
    for(int p2=0; p2<eps->get_len(); p2++){
        
        if(eps->get_value(p2)+eps->get_value(p1) <= eps->get_max_linspace()){
            p3_vals[p2] = new linspace_for_trap(0, eps->get_value(p2)+eps->get_value(p1), p2+p1+1);
            inner_vals[p2] = new dep_vars(p2+p1+1);
        }
        //p3 will just be eps
        else{
            p3_vals[p2] = new dummy_vars(eps);
            inner_vals[p2] = new dep_vars(eps->get_len());
        }
        
    }
    Fvv_values = new double**[4];
    Fvvbar_values = new double**[4];
    for(int i=0; i<4; i++){
        Fvv_values[i] = new double*[eps->get_len()];
        Fvvbar_values[i] = new double*[eps->get_len()];
        for(int j=0; j<eps->get_len(); j++){
            Fvv_values[i][j] = new double[eps->get_len()+1](); 
            Fvvbar_values[i][j] = new double[eps->get_len()+1]();  
        }
    }  
    
}

void nu_nu_collision_one::Fvvsc_components_term_1(density* dens, bool neutrino, int p2, int p3, double* F0, three_vector* F){
    
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();
    matrix* p_4 = new matrix();
    p_1->convert_p_to_identity_minus_matrix(dens, neutrino, p1);
    p_2->convert_p_to_identity_minus_matrix(dens, neutrino, p2);
    
    double max_lin = eps->get_max_linspace();
    double p3_energy = p3_vals[p2]->get_value(p3);
    //if p3 represents the last element in the p3_vals array and it is in the GL points, it must be interpolated
    p_3->convert_p_to_matrix(dens, neutrino, p3);
    
    
    double p4_energy = eps->get_value(p1)+eps->get_value(p2)-p3_energy;
    if(p4_energy<0){
        p4_energy = 0;
    }
    
    //this clause finds an interpolated value for the p4 matrix if p4_energy is not in the linspace
    if (eps->get_value(p1)<=max_lin and eps->get_value(p2)<=max_lin and p3_energy<=max_lin and p4_energy<=max_lin){
        p_4->convert_p_to_matrix(dens, neutrino, p1+p2-p3);
    }
    else{
        double A0 = dens->interpolate_p0(neutrino, p4_energy);
        three_vector* A = new three_vector();
        dens->interpolate_p0p(neutrino, p4_energy, A);
        p_4->convert_p_to_matrix(A0, A);
        delete A;
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



void nu_nu_collision_one::Fvvsc_components_term_2(density* dens, bool neutrino, int p2, int p3, double* F0, three_vector* F){
    
    
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();
    
    matrix* p_4 = new matrix(true);
    p_1->convert_p_to_matrix(dens, neutrino, p1);
    p_2->convert_p_to_matrix(dens, neutrino, p2);
    
    double max_lin = eps->get_max_linspace();
    double p3_energy = p3_vals[p2]->get_value(p3);
    
    
    //if p3 represents the last element in the p3_vals array and it is in the GL points, it must be interpolated
    p_3->convert_p_to_identity_minus_matrix(dens, neutrino, p3);
    
    
    double p4_energy = eps->get_value(p1)+eps->get_value(p2)-p3_energy;
    if(p4_energy<0){
        p4_energy = 0;
    }
    
    //this clause finds an interpolated value for the p4 matrix if p4_energy is not in the linspace
    if (eps->get_value(p1)<=max_lin and eps->get_value(p2)<=max_lin and p3_energy<=max_lin and p4_energy<=max_lin){
        p_4->convert_p_to_identity_minus_matrix(dens, neutrino, p1+p2-p3);
    }
    else{
        double A0 = dens->interpolate_p0(neutrino, p4_energy);
        three_vector* A = new three_vector();
        dens->interpolate_p0p(neutrino, p4_energy, A);
        p_4->convert_p_to_identity_minus_matrix(A0, A);
        delete A;
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

void nu_nu_collision_one::Fvvsc_components(density* dens, bool neutrino, int p2, int p3, double* F03, three_vector* F3){
    
    double F01;
    three_vector* F1 = new three_vector();
    double F02;
    three_vector* F2 = new three_vector();
    
    
    Fvvsc_components_term_1(dens, neutrino, p2, p3, &F01, F1);
    Fvvsc_components_term_2(dens, neutrino, p2, p3, &F02, F2);
    
    F2->multiply_by(-1);
    F3->add(F1, F2);
    
    *F03 = F01 - F02;
    
    delete F1;
    delete F2;
}

void nu_nu_collision_one::Fvvsc_for_p1(density* dens, bool neutrino){
    double F0 = 0;
    three_vector* Fxyz = new three_vector();
    for(int p2=0; p2<eps->get_len(); p2++){
        for(int p3=0; p3<p3_vals[p2]->get_len(); p3++){
            
            //only if p3_energy is less than p1_energy+p2_energy is this called
            if(p3_vals[p2]->get_value(p3) < eps->get_value(p1) + eps->get_value(p2)){
                if(eps->get_value(p1)+eps->get_value(p2)-p3_vals[p2]->get_value(p3)>=0){
                    Fvvsc_components(dens, neutrino, p2, p3, &F0, Fxyz);

                    Fvv_values[0][p2][p3] = F0;
                    Fvv_values[1][p2][p3] = Fxyz->get_value(0);
                    Fvv_values[2][p2][p3] = Fxyz->get_value(1);
                    Fvv_values[3][p2][p3] = Fxyz->get_value(2);
                }
            }
        }
    }  
    
    delete Fxyz;
}

void nu_nu_collision_one::Fvvbarsc_components_term_1(density* dens, bool neutrino, int p2, int p3, double* F0, three_vector* F){
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();
    
    matrix* p_4 = new matrix(true);
    p_1->convert_p_to_identity_minus_matrix(dens, neutrino, p1);
    p_2->convert_p_to_identity_minus_matrix(dens, not neutrino, p2);

    double max_lin = eps->get_max_linspace();
    double p3_energy = p3_vals[p2]->get_value(p3);
    
    p_3->convert_p_to_matrix(dens, neutrino, p3);
    
    double p4_energy = eps->get_value(p1)+eps->get_value(p2)-p3_energy;
    if(p4_energy<0){
        p4_energy = 0;
    }
    
    //this clause finds an interpolated value for the p4 matrix if p4_energy is not in the linspace
    if (eps->get_value(p1)<=max_lin and eps->get_value(p2)<=max_lin and eps->get_value(p3)<=max_lin and p4_energy<=max_lin){
        p_4->convert_p_to_matrix(dens, not neutrino, p1+p2-p3);
    }
    else{
        double A0 = dens->interpolate_p0(not neutrino, p4_energy);
        three_vector* A = new three_vector();
        dens->interpolate_p0p(not neutrino, p4_energy, A);
        p_4->convert_p_to_matrix(A0, A);
        delete A;
    }
    
    /*
    F_dummy1 = (rho_3)(1-rho_2)
    id1 = 1*tr((rho_3)(1-rho_2))
    F_dummy2 = (rho_3)(1-rho_2)+1*tr((rho_3)(1-rho_2))
    F_dummy3 = (1-rho_1)(rho_4)
    F_dummy4 = (1-rho_1)(rho_4) * [(rho_3)(1-rho_2)+1*tr((rho_3)(1-rho_2))]
    
    F_dummy5 = (rho_3)(rho_4)
    id2 = 1*tr((rho_3)(rho_4))
    F_dummy6 = (rho_3)(rho_4)+1*tr((rho_3)(rho_4))
    F_dummy7 = (1-rho_1)(1-rho_2)
    F_dummy8 = (1-rho_1)(1-rho_2) * [(rho_3)(rho_4)+1*tr((rho_3)(rho_4))]
    
    F_dummy9 = F_dummy4+F_dummy8
    */
    
    matrix* F_dummy1 = new matrix();
    F_dummy1->matrix_multiply(p_3, p_2);
    
    matrix* id1 = new matrix(true);
    id1->multiply_by(F_dummy1->get_A0()*(complex<double> (2,0)));
    
    matrix* F_dummy2 = new matrix();
    F_dummy2->matrix_add(F_dummy1, id1);
    
    matrix* F_dummy3 = new matrix();
    F_dummy3->matrix_multiply(p_1, p_4);
    
    matrix* F_dummy4 = new matrix();
    F_dummy4->matrix_multiply(F_dummy3, F_dummy2);
    
    matrix* F_dummy5 = new matrix();
    F_dummy5->matrix_multiply(p_3, p_4);
    
    matrix* id2 = new matrix(true);
    id2->multiply_by(F_dummy5->get_A0()*(complex<double> (2,0)));
    
    matrix* F_dummy6 = new matrix();
    F_dummy6->matrix_add(F_dummy5, id2);
    
    matrix* F_dummy7 = new matrix();
    F_dummy7->matrix_multiply(p_1, p_2);
    
    matrix* F_dummy8 = new matrix();
    F_dummy8->matrix_multiply(F_dummy7, F_dummy6);
    
    matrix* F_dummy9 = new matrix();
    F_dummy9->matrix_add(F_dummy4, F_dummy8);
    
    complex<double> comp_F0 = F_dummy9->get_A0();
    complex_three_vector* comp_F = F_dummy9->get_A();
    
    comp_F->multiply_by(2);
    
    *F0 = 2*real(comp_F0);
    F->make_real(comp_F);
    
    delete F_dummy1;
    delete F_dummy2;
    delete F_dummy3;
    delete F_dummy4;
    delete F_dummy5;
    delete F_dummy6;
    delete F_dummy7;
    delete F_dummy8;
    delete F_dummy9;
    delete id1;
    delete id2;
    delete p_1;
    delete p_2;
    delete p_3;
    delete p_4;    
}

void nu_nu_collision_one::Fvvbarsc_components_term_2(density* dens, bool neutrino, int p2, int p3, double* F0, three_vector* F){
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();
    
    matrix* p_4 = new matrix(true);
    p_1->convert_p_to_matrix(dens, neutrino, p1);
    p_2->convert_p_to_matrix(dens, not neutrino, p2);
    
    double max_lin = eps->get_max_linspace();
    double p3_energy = p3_vals[p2]->get_value(p3);
    
    p_3->convert_p_to_identity_minus_matrix(dens, neutrino, p3);
    
    double p4_energy = eps->get_value(p1)+eps->get_value(p2)-p3_energy;
    if(p4_energy<0){
        p4_energy = 0;
    }
    
    if (eps->get_value(p1)<=max_lin and eps->get_value(p2)<=max_lin and eps->get_value(p3)<=max_lin and p4_energy<=max_lin){
        p_4->convert_p_to_identity_minus_matrix(dens, not neutrino, p1+p2-p3);
    }
    else{
        double A0 = dens->interpolate_p0(not neutrino, p4_energy);
        three_vector* A = new three_vector();
        dens->interpolate_p0p(not neutrino, p4_energy, A);
        p_4->convert_p_to_identity_minus_matrix(A0, A);
        delete A;
    }
    
    
    /*
    F_dummy1 = (1-rho_3)(rho_2)
    id1 = 1*tr((1-rho_3)(rho_2))
    F_dummy2 = (1-rho_3)(rho_2)+1*tr((1-rho_3)(rho_2))
    F_dummy3 = (rho_1)(1-rho_4)
    F_dummy4 = (rho_1)(1-rho_4) * [(1-rho_3)(rho_2)+1*tr((1-rho_3)(rho_2))]
    
    F_dummy5 = (1-rho_3)(1-rho_4)
    id2 = 1*tr((1-rho_3)(1-rho_4))
    F_dummy6 = (1-rho_3)(1-rho_4)+1*tr((1-rho_3)(1-rho_4))
    F_dummy7 = (rho_1)(rho_2)
    F_dummy8 = (rho_1)(rho_2) * [(1-rho_3)(1-rho_4)+1*tr((1-rho_3)(1-rho_4))]
    
    F_dummy9 = F_dummy4+F_dummy8
    */
    
    matrix* F_dummy1 = new matrix();
    F_dummy1->matrix_multiply(p_3, p_2);
    
    matrix* id1 = new matrix(true);
    id1->multiply_by(F_dummy1->get_A0()*(complex<double> (2,0)));
    
    matrix* F_dummy2 = new matrix();
    F_dummy2->matrix_add(F_dummy1, id1);
    
    matrix* F_dummy3 = new matrix();
    F_dummy3->matrix_multiply(p_1, p_4);
    
    matrix* F_dummy4 = new matrix();
    F_dummy4->matrix_multiply(F_dummy3, F_dummy2);
    
    matrix* F_dummy5 = new matrix();
    F_dummy5->matrix_multiply(p_3, p_4);
    
    matrix* id2 = new matrix(true);
    id2->multiply_by(F_dummy5->get_A0()*(complex<double> (2,0)));
    
    matrix* F_dummy6 = new matrix();
    F_dummy6->matrix_add(F_dummy5, id2);
    
    matrix* F_dummy7 = new matrix();
    F_dummy7->matrix_multiply(p_1, p_2);
    
    matrix* F_dummy8 = new matrix();
    F_dummy8->matrix_multiply(F_dummy7, F_dummy6);
    
    matrix* F_dummy9 = new matrix();
    F_dummy9->matrix_add(F_dummy4, F_dummy8);
    
    complex<double> comp_F0 = F_dummy9->get_A0();
    complex_three_vector* comp_F = F_dummy9->get_A();
    
    comp_F->multiply_by(2);
    
    *F0 = 2*real(comp_F0);
    F->make_real(comp_F);
    
    delete F_dummy1;
    delete F_dummy2;
    delete F_dummy3;
    delete F_dummy4;
    delete F_dummy5;
    delete F_dummy6;
    delete F_dummy7;
    delete F_dummy8;
    delete F_dummy9;
    delete id1;
    delete id2;
    delete p_1;
    delete p_2;
    delete p_3;
    delete p_4;
}

void nu_nu_collision_one::Fvvbarsc_components(density* dens, bool neutrino, int p2, int p3, double* F03, three_vector* F3){
       
    double F01;
    three_vector* F1 = new three_vector();
    double F02;
    three_vector* F2 = new three_vector();
    
    Fvvbarsc_components_term_1(dens, neutrino, p2, p3, &F01, F1);
    Fvvbarsc_components_term_2(dens, neutrino, p2, p3, &F02, F2);
    
    F2->multiply_by(-1);
    F3->add(F1, F2);
    
    *F03 = F01 - F02;
    
    delete F1;
    delete F2;
}

void nu_nu_collision_one::Fvvbarsc_for_p1(density* dens, bool neutrino){
    double F0 = 0;
    three_vector* Fxyz = new three_vector();
    for(int p2=0; p2<eps->get_len(); p2++){
        for(int p3=0; p3<p3_vals[p2]->get_len(); p3++){
            
            //only if p3_energy is less than p1_energy+p2_energy is this called--> will make integrand 0 past p1+p2
            if(p3_vals[p2]->get_value(p3) < eps->get_value(p1) + eps->get_value(p2)){
                if(eps->get_value(p1)+eps->get_value(p2)-p3_vals[p2]->get_value(p3)>=0){
                    Fvvbarsc_components(dens, neutrino, p2, p3, &F0, Fxyz);

                    Fvvbar_values[0][p2][p3] = F0;
                    Fvvbar_values[1][p2][p3] = Fxyz->get_value(0);
                    Fvvbar_values[2][p2][p3] = Fxyz->get_value(1);
                    Fvvbar_values[3][p2][p3] = Fxyz->get_value(2);
                    
                }
            }
        }
    }  
    
    delete Fxyz;
}


double nu_nu_collision_one::J1(double p1, double p2, double p3){
    return 16./15 * pow(p3,3) * (10 * pow(p1+p2,2) - 15 * (p1+p2) * p3 + 6*pow(p3,2));  
}

double nu_nu_collision_one::J2(double p1, double p2){
    return 16./15 * pow(p2,3) * (10 * pow(p1,2) + 5 * p1*p2 + pow(p2,2));  
}

double nu_nu_collision_one::J3(double p1, double p2, double p3){
    return 16./15 * (pow(p1+p2,5) - 10 * pow(p1+p2,2) * pow(p3, 3) + 15 * (p1+p2) * pow(p3,4) - 6 * pow(p3,5));
}

double nu_nu_collision_one::K1(double p1, double p3){
    return 16./15 * pow(p3,3) * (10 * pow(p1,2) - 5 * p1*p3 + pow(p3,2));
}

double nu_nu_collision_one::K2(double p1, double p2, double p3){
    return 16./15 * pow(p2,3) * (10 * pow(p1-p3,2) + 15 * (p1-p3) * p2 + 6 * pow(p2,2));
}

double nu_nu_collision_one::K3(double p1, double p2, double p3){
    return 16./15 * (pow(p1-p3,5) + 10 * pow(p1-p3,2) * pow(p2,3) + 15 * (p1-p3) * pow(p2,4) + 6 * pow(p2,5));
}


double nu_nu_collision_one::interior_integral(int p2, int which_term){
    double p_1_energy = eps->get_value(p1);
    double max_energy = eps->get_value(eps->get_len()-1);
    
    if(p2<p1){
        for(int p3=0; p3<p2; p3++){
            inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J1(p_1_energy, eps->get_value(p2), p3_vals[p2]->get_value(p3)) + Fvvbar_values[which_term][p2][p3] * K1(p_1_energy, p3_vals[p2]->get_value(p3)));
        }
        for(int p3=p2; p3<p1; p3++){
            inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J2(p_1_energy, eps->get_value(p2)) + Fvvbar_values[which_term][p2][p3] * K2(p_1_energy, eps->get_value(p2), p3_vals[p2]->get_value(p3)));
        }
        for(int p3=p1; p3<p3_vals[p2]->get_len(); p3++){
            inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J3(p_1_energy, eps->get_value(p2), p3_vals[p2]->get_value(p3)) + Fvvbar_values[which_term][p2][p3] * K3(p_1_energy, eps->get_value(p2), p3_vals[p2]->get_value(p3)));
        }

    }

    else{
        for(int p3=0; p3<p1; p3++){
            inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J1(p_1_energy, eps->get_value(p2), p3_vals[p2]->get_value(p3)) + Fvvbar_values[which_term][p2][p3] * K1(p_1_energy, p3_vals[p2]->get_value(p3)));
        }
        for(int p3=p1; p3<p2; p3++){
            inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J2(eps->get_value(p2), p_1_energy) + Fvvbar_values[which_term][p2][p3] * K1(p3_vals[p2]->get_value(p3), p_1_energy));

        }
        for(int p3=p2; p3<p3_vals[p2]->get_len(); p3++){
            inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J3(p_1_energy, eps->get_value(p2), p3_vals[p2]->get_value(p3)) + Fvvbar_values[which_term][p2][p3] * K3(p_1_energy, eps->get_value(p2), p3_vals[p2]->get_value(p3)));
        }

    }
    
    if(p2==200){
        std::cout << "p2=" << p2 << std::endl;
        for(int p3=0; p3<p3_vals[p2]->get_len(); p3++){
            std::cout << Fvvbar_values[0][p2][p3] << ", ";
        }
        std::cout << std::endl << "------" << std::endl;
        
        for(int i=0; i<p3_vals[p2]->get_len(); i++){
            std::cout << p3_vals[p2]->get_value(i) << ", ";
        }
    }
    
    double result = p3_vals[p2]->integrate(inner_vals[p2]);
    
    return result;
}

//note: results must be length 4
void nu_nu_collision_one::whole_integral(density* dens, bool neutrino, double* results){
    if (p1==0){
        for(int i=0; i<4; i++){
            results[i] = 0;
        }
    }
    else{
        //populates F_values
        Fvvsc_for_p1(dens, neutrino);
        Fvvbarsc_for_p1(dens, neutrino);
        double Tcm = dens->get_Tcm();
            
        double p_1_energy = eps->get_value(p1);
        for(int i=0; i<4; i++){
            for(int p2=0; p2<eps->get_len(); p2++){
                outer_vals->set_value(p2, interior_integral(p2, i));
            }
            results[i] = eps->integrate(outer_vals);
            results[i] *= pow(Tcm, 5) * pow(_GF_,2) / (pow(2*_PI_,3) * pow(p_1_energy,2));
        }
    }
}


nu_nu_collision_one::~nu_nu_collision_one(){
    delete outer_vals;
    for(int i=0; i<eps->get_len(); i++){
        delete inner_vals[i];
        delete p3_vals[i];
    }
    delete[] inner_vals;
    delete[] p3_vals;
    for(int i=0; i<4; i++){
        for(int j=0; j<eps->get_len(); j++){
            delete[] Fvv_values[i][j];
            delete[] Fvvbar_values[i][j];
        }
        delete[] Fvv_values[i];
        delete[] Fvvbar_values[i];
    }
    delete[] Fvv_values;
    delete[] Fvvbar_values;
    delete eps;
    
}


nu_nu_collision_two::nu_nu_collision_two(linspace_and_gl* e, int p1_index){
    eps = new linspace_and_gl(e);
    p1 = p1_index;
    count = 0;
    outer_vals = new dep_vars(eps->get_len());
    inner_vals = new dep_vars*[eps->get_len()];
    p3_vals = new dummy_vars*[eps->get_len()];
    for(int p2=0; p2<eps->get_len(); p2++){
        
        if(eps->get_value(p2)+eps->get_value(p1) <= eps->get_max_linspace()){
            p3_vals[p2] = new linspace_for_trap(0, eps->get_value(p2)+eps->get_value(p1), p2+p1+1);
            inner_vals[p2] = new dep_vars(p2+p1+1);
        }
        else{
            p3_vals[p2] = new linspace_and_gel(eps, eps->get_value(p2)+eps->get_value(p1), 10);
            
            inner_vals[p2] = new dep_vars(p3_vals[p2]->get_len());
        }
    }
    Fvv_values = new double**[4];
    Fvvbar_values = new double**[4];
    for(int i=0; i<4; i++){
        Fvv_values[i] = new double*[eps->get_len()];
        Fvvbar_values[i] = new double*[eps->get_len()];
        for(int j=0; j<eps->get_len(); j++){
            Fvv_values[i][j] = new double[eps->get_len()+1](); 
            Fvvbar_values[i][j] = new double[eps->get_len()+1]();  
        }
    }  
}


void nu_nu_collision_two::Fvvsc_components_term_1(density* dens, bool neutrino, int p2, int p3, double* F0, three_vector* F){
    
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();
    matrix* p_4 = new matrix();
    p_1->convert_p_to_identity_minus_matrix(dens, neutrino, p1);
    p_2->convert_p_to_identity_minus_matrix(dens, neutrino, p2);
    
    double max_lin = eps->get_max_linspace();
    double p3_energy = p3_vals[p2]->get_value(p3);
    
    //if p3 represents the last element in the p3_vals array and it is in the GL points, it must be interpolated
    if(p3_energy>max_lin){
        double A0 = dens->interpolate_p0(neutrino, p3_energy);
        three_vector* A = new three_vector();
        dens->interpolate_p0p(neutrino, p3_energy, A);
        p_3->convert_p_to_matrix(A0, A);
        delete A;
    }
    else{
        p_3->convert_p_to_matrix(dens, neutrino, p3);
    }
    
    double p4_energy = eps->get_value(p1)+eps->get_value(p2)-p3_energy;
    //this prevents computer errors in subtraction; if p3 is p1+p2 p4 stays identically zero
    if(p4_energy<0){
        p4_energy = 0;
    }
    //this clause finds an interpolated value for the p4 matrix if p4_energy is not in the linspace
    if (eps->get_value(p1)<=max_lin and eps->get_value(p2)<=max_lin and eps->get_value(p3)<=max_lin and p4_energy<=max_lin){
            p_4->convert_p_to_matrix(dens, neutrino, p1+p2-p3);
    }
    else{
        double A0 = dens->interpolate_p0(neutrino, p4_energy);
        three_vector* A = new three_vector();
        dens->interpolate_p0p(neutrino, p4_energy, A);
        p_4->convert_p_to_matrix(A0, A);
        delete A;
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


void nu_nu_collision_two::Fvvsc_components_term_2(density* dens, bool neutrino, int p2, int p3, double* F0, three_vector* F){
    
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();
    
    matrix* p_4 = new matrix(true);
    p_1->convert_p_to_matrix(dens, neutrino, p1);
    p_2->convert_p_to_matrix(dens, neutrino, p2);
    
    double max_lin = eps->get_max_linspace();
    double p3_energy = p3_vals[p2]->get_value(p3);
    
    
    //if p3 represents the last element in the p3_vals array and it is in the GL points, it must be interpolated
    if(p3_energy>max_lin){
        double A0 = dens->interpolate_p0(neutrino, p3_energy);
        three_vector* A = new three_vector();
        dens->interpolate_p0p(neutrino, p3_energy, A);
        p_3->convert_p_to_identity_minus_matrix(A0, A);
        delete A;
    }
    else{
        p_3->convert_p_to_identity_minus_matrix(dens, neutrino, p3);
    }
    
    double p4_energy = eps->get_value(p1)+eps->get_value(p2)-p3_energy;
    //this prevents computer errors in subtraction; if p3 is p1+p2 p4 stays identically zero
    if(p4_energy<0){
        p4_energy = 0;
    }
    //this clause finds an interpolated value for the p4 matrix if p4_energy is not in the linspace
    if (eps->get_value(p1)<=max_lin and eps->get_value(p2)<=max_lin and eps->get_value(p3)<=max_lin and p4_energy<=max_lin){
        p_4->convert_p_to_identity_minus_matrix(dens, neutrino, p1+p2-p3);
    }
    else{
        double A0 = dens->interpolate_p0(neutrino, p4_energy);
        three_vector* A = new three_vector();
        dens->interpolate_p0p(neutrino, p4_energy, A);
        p_4->convert_p_to_identity_minus_matrix(A0, A);
        delete A;
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

void nu_nu_collision_two::Fvvsc_components(density* dens, bool neutrino, int p2, int p3, double* F03, three_vector* F3){
    
    double F01;
    three_vector* F1 = new three_vector();
    double F02;
    three_vector* F2 = new three_vector();
    
    
    Fvvsc_components_term_1(dens, neutrino, p2, p3, &F01, F1);
    Fvvsc_components_term_2(dens, neutrino, p2, p3, &F02, F2);
    
    F2->multiply_by(-1);
    F3->add(F1, F2);
    
    *F03 = F01 - F02;
    
    delete F1;
    delete F2;
}

void nu_nu_collision_two::Fvvsc_for_p1(density* dens, bool neutrino){
    double F0 = 0;
    three_vector* Fxyz = new three_vector();
    for(int p2=0; p2<eps->get_len(); p2++){
        for(int p3=0; p3<p3_vals[p2]->get_len(); p3++){
            
            if(eps->get_value(p1)+eps->get_value(p2)-p3_vals[p2]->get_value(p3)>=0){                
                Fvvsc_components(dens, neutrino, p2, p3, &F0, Fxyz);
                
                Fvv_values[0][p2][p3] = F0;
                Fvv_values[1][p2][p3] = Fxyz->get_value(0);
                Fvv_values[2][p2][p3] = Fxyz->get_value(1);
                Fvv_values[3][p2][p3] = Fxyz->get_value(2);
            }
        }
    }  
    
    delete Fxyz;
}

void nu_nu_collision_two::Fvvbarsc_components_term_1(density* dens, bool neutrino, int p2, int p3, double* F0, three_vector* F){
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();
    
    matrix* p_4 = new matrix(true);
    p_1->convert_p_to_identity_minus_matrix(dens, neutrino, p1);
    p_2->convert_p_to_identity_minus_matrix(dens, not neutrino, p2);

    double max_lin = eps->get_max_linspace();
    double p3_energy = p3_vals[p2]->get_value(p3);
    
    
    //if p3 represents the last element in the p3_vals array and it is in the GL points, it must be interpolated
    if(p3_energy>max_lin){
        double A0 = dens->interpolate_p0(neutrino, p3_energy);
        three_vector* A = new three_vector();
        dens->interpolate_p0p(neutrino, p3_energy, A);
        p_3->convert_p_to_matrix(A0, A);
        delete A;
    }
    else{
        p_3->convert_p_to_matrix(dens, neutrino, p3);
    }
    
    double p4_energy = eps->get_value(p1)+eps->get_value(p2)-p3_energy;
    //this prevents computer errors in subtraction; if p3 is p1+p2 p4 stays identically zero
    if(p4_energy<0){
        p4_energy = 0;
    }
    
    //this clause finds an interpolated value for the p4 matrix if p4_energy is not in the linspace
    if (eps->get_value(p1)<=max_lin and eps->get_value(p2)<=max_lin and eps->get_value(p3)<=max_lin and p4_energy<=max_lin){
        p_4->convert_p_to_matrix(dens, not neutrino, p1+p2-p3);
    }
    else{
        double A0 = dens->interpolate_p0(not neutrino, p4_energy);
        three_vector* A = new three_vector();
        dens->interpolate_p0p(not neutrino, p4_energy, A);
        p_4->convert_p_to_matrix(A0, A);
        delete A;
    }
    
    /*
    F_dummy1 = (rho_3)(1-rho_2)
    id1 = 1*tr((rho_3)(1-rho_2))
    F_dummy2 = (rho_3)(1-rho_2)+1*tr((rho_3)(1-rho_2))
    F_dummy3 = (1-rho_1)(rho_4)
    F_dummy4 = (1-rho_1)(rho_4) * [(rho_3)(1-rho_2)+1*tr((rho_3)(1-rho_2))]
    
    F_dummy5 = (rho_3)(rho_4)
    id2 = 1*tr((rho_3)(rho_4))
    F_dummy6 = (rho_3)(rho_4)+1*tr((rho_3)(rho_4))
    F_dummy7 = (1-rho_1)(1-rho_2)
    F_dummy8 = (1-rho_1)(1-rho_2) * [(rho_3)(rho_4)+1*tr((rho_3)(rho_4))]
    
    F_dummy9 = F_dummy4+F_dummy8
    */
    
    matrix* F_dummy1 = new matrix();
    F_dummy1->matrix_multiply(p_3, p_2);
    
    matrix* id1 = new matrix(true);
    id1->multiply_by(F_dummy1->get_A0()*(complex<double> (2,0)));
    
    matrix* F_dummy2 = new matrix();
    F_dummy2->matrix_add(F_dummy1, id1);
    
    matrix* F_dummy3 = new matrix();
    F_dummy3->matrix_multiply(p_1, p_4);
    
    matrix* F_dummy4 = new matrix();
    F_dummy4->matrix_multiply(F_dummy3, F_dummy2);
    
    matrix* F_dummy5 = new matrix();
    F_dummy5->matrix_multiply(p_3, p_4);
    
    matrix* id2 = new matrix(true);
    id2->multiply_by(F_dummy5->get_A0()*(complex<double> (2,0)));
    
    matrix* F_dummy6 = new matrix();
    F_dummy6->matrix_add(F_dummy5, id2);
    
    matrix* F_dummy7 = new matrix();
    F_dummy7->matrix_multiply(p_1, p_2);
    
    matrix* F_dummy8 = new matrix();
    F_dummy8->matrix_multiply(F_dummy7, F_dummy6);
    
    matrix* F_dummy9 = new matrix();
    F_dummy9->matrix_add(F_dummy4, F_dummy8);
    
    complex<double> comp_F0 = F_dummy9->get_A0();
    complex_three_vector* comp_F = F_dummy9->get_A();
    
    comp_F->multiply_by(2);
    
    *F0 = 2*real(comp_F0);
    F->make_real(comp_F);
    
    delete F_dummy1;
    delete F_dummy2;
    delete F_dummy3;
    delete F_dummy4;
    delete F_dummy5;
    delete F_dummy6;
    delete F_dummy7;
    delete F_dummy8;
    delete F_dummy9;
    delete id1;
    delete id2;
    delete p_1;
    delete p_2;
    delete p_3;
    delete p_4;    
}

void nu_nu_collision_two::Fvvbarsc_components_term_2(density* dens, bool neutrino, int p2, int p3, double* F0, three_vector* F){
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();
    
    matrix* p_4 = new matrix(true);
    p_1->convert_p_to_matrix(dens, neutrino, p1);
    p_2->convert_p_to_matrix(dens, not neutrino, p2);
    
    double max_lin = eps->get_max_linspace();
    double p3_energy = p3_vals[p2]->get_value(p3);
    
    //if p3 represents the last element in the p3_vals array and it is in the GL points, it must be interpolated
    if(p3_energy>max_lin){
        double A0 = dens->interpolate_p0(neutrino, p3_energy);
        three_vector* A = new three_vector();
        dens->interpolate_p0p(neutrino, p3_energy, A);
        p_3->convert_p_to_identity_minus_matrix(A0, A);
        delete A;
        
    }
    else{
        p_3->convert_p_to_identity_minus_matrix(dens, neutrino, p3);
    }
    
    double p4_energy = eps->get_value(p1)+eps->get_value(p2)-p3_energy;
    //this prevents computer errors in subtraction; if p3 is p1+p2 p4 stays identically zero
    if(p4_energy<0){
        p4_energy = 0;
    }
    
    //this clause finds an interpolated value for the p4 matrix if p4_energy is not in the linspace
    if (eps->get_value(p1)<=max_lin and eps->get_value(p2)<=max_lin and eps->get_value(p3)<=max_lin and p4_energy<=max_lin){
        p_4->convert_p_to_identity_minus_matrix(dens, not neutrino, p1+p2-p3);
    }
    else{
        double A0 = dens->interpolate_p0(not neutrino, p4_energy);
        three_vector* A = new three_vector();
        dens->interpolate_p0p(not neutrino, p4_energy, A);
        p_4->convert_p_to_identity_minus_matrix(A0, A);
        delete A;
    }
    
    /*
    F_dummy1 = (1-rho_3)(rho_2)
    id1 = 1*tr((1-rho_3)(rho_2))
    F_dummy2 = (1-rho_3)(rho_2)+1*tr((1-rho_3)(rho_2))
    F_dummy3 = (rho_1)(1-rho_4)
    F_dummy4 = (rho_1)(1-rho_4) * [(1-rho_3)(rho_2)+1*tr((1-rho_3)(rho_2))]
    
    F_dummy5 = (1-rho_3)(1-rho_4)
    id2 = 1*tr((1-rho_3)(1-rho_4))
    F_dummy6 = (1-rho_3)(1-rho_4)+1*tr((1-rho_3)(1-rho_4))
    F_dummy7 = (rho_1)(rho_2)
    F_dummy8 = (rho_1)(rho_2) * [(1-rho_3)(1-rho_4)+1*tr((1-rho_3)(1-rho_4))]
    
    F_dummy9 = F_dummy4+F_dummy8
    */
    
    matrix* F_dummy1 = new matrix();
    F_dummy1->matrix_multiply(p_3, p_2);
    
    matrix* id1 = new matrix(true);
    id1->multiply_by(F_dummy1->get_A0()*(complex<double> (2,0)));
    
    matrix* F_dummy2 = new matrix();
    F_dummy2->matrix_add(F_dummy1, id1);
    
    matrix* F_dummy3 = new matrix();
    F_dummy3->matrix_multiply(p_1, p_4);
    
    matrix* F_dummy4 = new matrix();
    F_dummy4->matrix_multiply(F_dummy3, F_dummy2);
    
    matrix* F_dummy5 = new matrix();
    F_dummy5->matrix_multiply(p_3, p_4);
    
    matrix* id2 = new matrix(true);
    id2->multiply_by(F_dummy5->get_A0()*(complex<double> (2,0)));
    
    matrix* F_dummy6 = new matrix();
    F_dummy6->matrix_add(F_dummy5, id2);
    
    matrix* F_dummy7 = new matrix();
    F_dummy7->matrix_multiply(p_1, p_2);
    
    matrix* F_dummy8 = new matrix();
    F_dummy8->matrix_multiply(F_dummy7, F_dummy6);
    
    matrix* F_dummy9 = new matrix();
    F_dummy9->matrix_add(F_dummy4, F_dummy8);
    
    complex<double> comp_F0 = F_dummy9->get_A0();
    complex_three_vector* comp_F = F_dummy9->get_A();
    
    comp_F->multiply_by(2);
    
    *F0 = 2*real(comp_F0);
    F->make_real(comp_F);
    
    delete F_dummy1;
    delete F_dummy2;
    delete F_dummy3;
    delete F_dummy4;
    delete F_dummy5;
    delete F_dummy6;
    delete F_dummy7;
    delete F_dummy8;
    delete F_dummy9;
    delete id1;
    delete id2;
    delete p_1;
    delete p_2;
    delete p_3;
    delete p_4;
}

void nu_nu_collision_two::Fvvbarsc_components(density* dens, bool neutrino, int p2, int p3, double* F03, three_vector* F3){
    
    double F01;
    three_vector* F1 = new three_vector();
    double F02;
    three_vector* F2 = new three_vector();
    
    
    Fvvbarsc_components_term_1(dens, neutrino, p2, p3, &F01, F1);
    Fvvbarsc_components_term_2(dens, neutrino, p2, p3, &F02, F2);
    
    F2->multiply_by(-1);
    F3->add(F1, F2);
    
    *F03 = F01 - F02;
    
    delete F1;
    delete F2;
}

void nu_nu_collision_two::Fvvbarsc_for_p1(density* dens, bool neutrino){
    double F0 = 0;
    three_vector* Fxyz = new three_vector();
    for(int p2=0; p2<eps->get_len(); p2++){
        for(int p3=0; p3<p3_vals[p2]->get_len(); p3++){
            
            if(eps->get_value(p1)+eps->get_value(p2)-p3_vals[p2]->get_value(p3)>=0){
                Fvvbarsc_components(dens, neutrino, p2, p3, &F0, Fxyz);
                
                Fvvbar_values[0][p2][p3] = F0;
                Fvvbar_values[1][p2][p3] = Fxyz->get_value(0);
                Fvvbar_values[2][p2][p3] = Fxyz->get_value(1);
                Fvvbar_values[3][p2][p3] = Fxyz->get_value(2);
            }
        }
    }  
    
    delete Fxyz;
}


double nu_nu_collision_two::J1(double p1, double p2, double p3){
    return 16./15 * pow(p3,3) * (10 * pow(p1+p2,2) - 15 * (p1+p2) * p3 + 6*pow(p3,2));  
}

double nu_nu_collision_two::J2(double p1, double p2){
    return 16./15 * pow(p2,3) * (10 * pow(p1,2) + 5 * p1*p2 + pow(p2,2));  
}

double nu_nu_collision_two::J3(double p1, double p2, double p3){
    return 16./15 * (pow(p1+p2,5) - 10 * pow(p1+p2,2) * pow(p3, 3) + 15 * (p1+p2) * pow(p3,4) - 6 * pow(p3,5));
}

double nu_nu_collision_two::K1(double p1, double p3){
    return 16./15 * pow(p3,3) * (10 * pow(p1,2) - 5 * p1*p3 + pow(p3,2));
}

double nu_nu_collision_two::K2(double p1, double p2, double p3){
    return 16./15 * pow(p2,3) * (10 * pow(p1-p3,2) + 15 * (p1-p3) * p2 + 6 * pow(p2,2));
}

double nu_nu_collision_two::K3(double p1, double p2, double p3){
    return 16./15 * (pow(p1-p3,5) + 10 * pow(p1-p3,2) * pow(p2,3) + 15 * (p1-p3) * pow(p2,4) + 6 * pow(p2,5));
}


double nu_nu_collision_two::interior_integral(int p2, int which_term){
    double p_1_energy = eps->get_value(p1);
    double max_energy = eps->get_value(eps->get_len()-1);
    

    if(p2<p1){
        for(int p3=0; p3<p2; p3++){
            inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J1(p_1_energy, eps->get_value(p2), p3_vals[p2]->get_value(p3)) + Fvvbar_values[which_term][p2][p3] * K1(p_1_energy, p3_vals[p2]->get_value(p3)));
        }
        for(int p3=p2; p3<p1; p3++){
            inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J2(p_1_energy, eps->get_value(p2)) + Fvvbar_values[which_term][p2][p3] * K2(p_1_energy, eps->get_value(p2), p3_vals[p2]->get_value(p3)));
        }
        for(int p3=p1; p3<p3_vals[p2]->get_len(); p3++){
            inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J3(p_1_energy, eps->get_value(p2), p3_vals[p2]->get_value(p3)) + Fvvbar_values[which_term][p2][p3] * K3(p_1_energy, eps->get_value(p2), p3_vals[p2]->get_value(p3)));
        }
        
    }
    
    else{
        for(int p3=0; p3<p1; p3++){
            inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J1(p_1_energy, eps->get_value(p2), p3_vals[p2]->get_value(p3)) + Fvvbar_values[which_term][p2][p3] * K1(p_1_energy, p3_vals[p2]->get_value(p3)));
        }
        for(int p3=p1; p3<p2; p3++){
            inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J2(eps->get_value(p2), p_1_energy) + Fvvbar_values[which_term][p2][p3] * K1(p3_vals[p2]->get_value(p3), p_1_energy));
            
        }
        for(int p3=p2; p3<p3_vals[p2]->get_len(); p3++){
            inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J3(p_1_energy, eps->get_value(p2), p3_vals[p2]->get_value(p3)) + Fvvbar_values[which_term][p2][p3] * K3(p_1_energy, eps->get_value(p2), p3_vals[p2]->get_value(p3)));
        }
        
    }
    /*
    if(p2==150){
        std::cout << "p2=" << p2 << std::endl;
        for(int i=0; i<inner_vals[p2]->length(); i++){
            std::cout << inner_vals[p2]->get_value(i) << ", ";
        }
        std::cout << std::endl << "------" << std::endl;
        
        for(int i=0; i<p3_vals[p2]->get_len(); i++){
            std::cout << p3_vals[p2]->get_value(i) << ", ";
        }
    }*/
    double result = p3_vals[p2]->integrate(inner_vals[p2]);
    
    
    return result;
}

//note: results must be length 4
void nu_nu_collision_two::whole_integral(density* dens, bool neutrino, double* results){
    if (p1==0){
        for(int i=0; i<4; i++){
            results[i] = 0;
        }
    }
    else{
        //populates F_values
        Fvvsc_for_p1(dens, neutrino);
        Fvvbarsc_for_p1(dens, neutrino);
        double Tcm = dens->get_Tcm();
            
        double p_1_energy = eps->get_value(p1);
        for(int i=0; i<4; i++){
            for(int p2=0; p2<eps->get_len(); p2++){
                outer_vals->set_value(p2, interior_integral(p2, i));
            }
            results[i] = eps->integrate(outer_vals);
            results[i] *= pow(Tcm, 5) * pow(_GF_,2) / (pow(2*_PI_,3) * pow(p_1_energy,2));
        }
    }
}


nu_nu_collision_two::~nu_nu_collision_two(){
    delete outer_vals;
    for(int i=0; i<eps->get_len(); i++){
        delete inner_vals[i];
        delete p3_vals[i];
    }
    delete[] inner_vals;
    delete[] p3_vals;
    for(int i=0; i<4; i++){
        for(int j=0; j<eps->get_len(); j++){
            delete[] Fvv_values[i][j];
            delete[] Fvvbar_values[i][j];
        }
        delete[] Fvv_values[i];
        delete[] Fvvbar_values[i];
    }
    delete[] Fvv_values;
    delete[] Fvvbar_values;
    delete eps;
    
}
