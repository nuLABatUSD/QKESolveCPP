#include "QKE_methods.hh"
#include <iostream>
#include "constants.hh"
#include <cmath>
#include "thermodynamics.hh"
#include "gl_vals.hh"
#include "matrices.hh"


void three_vector_for_QKE::v_vacuum(double delta_m_squared, double cos_2theta, double sin_2theta ){
    values[0] = delta_m_squared / 2. * sin_2theta;
    values[1] = 0.;
    values[2] = - delta_m_squared / 2. * cos_2theta;
}


void three_vector_for_QKE::v_thermal(dummy_vars* q, density* d){
    
    dep_vars* d0 = new dep_vars(q->get_len()); 
    dep_vars* d1 = new dep_vars(q->get_len()); 
    dep_vars* d2 = new dep_vars(q->get_len());
    three_vector* dummy1 = new three_vector();
    three_vector* dummy2 = new three_vector();
    
    for(int i=0; i<q->get_len(); i++){
        d->p0_p(i, true, dummy1);
        d->p0_p(i, false, dummy2);

        d0->set_value(i, pow(q->get_value(i),3) * (dummy1->get_value(0) + dummy2->get_value(0)));
        d1->set_value(i, pow(q->get_value(i),3) * (dummy1->get_value(1) + dummy2->get_value(1)));
        d2->set_value(i, pow(q->get_value(i),3) * (dummy1->get_value(2) + dummy2->get_value(2))); 
    }
    values[0] = q->integrate(d0);
    values[1] = q->integrate(d1);
    values[2] = q->integrate(d2);
                      
    delete d0;
    delete d1;
    delete d2;
    delete dummy1;
    delete dummy2;
   
    
    for(int i=0; i<3; i++){
        values[i] *= -8*sqrt(2)*_GF_/(3*pow(_Z_boson_,2)) * pow(d->get_Tcm(), 4);
    }
    
    double energy_dens = 0;
    double pressure = 0;
    energy_and_pressure(_electron_mass_, d->get_T(), &energy_dens, &pressure);
    values[2] += -2 * sqrt(2) * pow(_W_boson_,-2) * _GF_ * (energy_dens+pressure);
    
}

void three_vector_for_QKE::v_density(dummy_vars* q, density* d){
    
    
    dep_vars* d0 = new dep_vars(q->get_len()); 
    dep_vars* d1 = new dep_vars(q->get_len()); 
    dep_vars* d2 = new dep_vars(q->get_len()); 
    three_vector* dummy1 = new three_vector();
    three_vector* dummy2 = new three_vector();

    for (int i=0; i<q->get_len(); i++){
        d->p0_p(i, true, dummy1);
        d->p0_p(i, false, dummy2);
        d0->set_value(i,pow(q->get_value(i),2) * (dummy1->get_value(0) - dummy2->get_value(0)));
        d1->set_value(i,pow(q->get_value(i),2) * (dummy1->get_value(1) - dummy2->get_value(1)));
        d2->set_value(i,pow(q->get_value(i),2) * (dummy1->get_value(2) - dummy2->get_value(2)));
    }
    values[0] = q->integrate(d0);
    values[1] = q->integrate(d1);
    values[2] = q->integrate(d2);
                      
    delete d0;
    delete d1;
    delete d2;
    delete dummy1;
    delete dummy2;

    
    for (int i=0; i<3; i++){
        values[i] *= sqrt(2)*_GF_ / (2 * pow(_PI_,2)) * pow(d->get_Tcm(), 3);
    }
}

//density
density::density(int num, dummy_vars* eps):dep_vars(8*num+2)
{
    N_bins = num;
    E = new dummy_vars(eps);
    
}

density::density(int num, dummy_vars* eps, double* dvals):dep_vars(8*num+2){
    N_bins = num;
    E = new dummy_vars(eps);
    for(int i=0; i<N_bins*8+2; i++){
        values[i] = dvals[i];
    }
    
    
}

density::density(dummy_vars* eps, double eta_nu, double eta_mu):dep_vars(8*eps->get_len()+2)
{
    N_bins = eps->get_len();
    E = new dummy_vars(eps);

    double fnu = 0;
    double fmu = 0;
    double fnubar = 0;
    double fmubar = 0;

    double eps_temp = 0.;
    
    for (int i=0; i<N_bins; i++){
        eps_temp = eps->get_value(i);
        fnu = 1 / (exp(eps_temp - eta_nu)+1);
        fmu = 1 / (exp(eps_temp - eta_mu)+1);
        values[4*i] = fnu + fmu;
        values[4*i+3] =  (fnu - fmu)/(fnu+fmu);
       
       fnubar = 1 / (exp(eps_temp + eta_nu)+1);
       fmubar = 1 / (exp(eps_temp + eta_mu)+1);
       values[4*N_bins + 4*i] = fnubar + fmubar;
       values[4*N_bins + 4*i+3] = (fnubar - fmubar)/(fnubar+fmubar);
      
    }
    
}

density::density(density* copy_me):dep_vars(copy_me)
{
    
    N_bins = copy_me->num_bins();
    E = new dummy_vars(copy_me->get_E());

}

density::~density()
{    delete E; }

dummy_vars* density::get_E(){
    return E;
}

double density::get_E_value(int i){
    return E->get_value(i);
}

double density::get_T(){
    return values[N-2];
}

double density::get_Tcm()
{return values[N-1];}

int density::num_bins()
{return N_bins;}

void density::set_T(double T)
{ 
    values[N-2] = T;
    values[N-1] = T;
}

double density::p0(int t, bool neutrino){
    if(neutrino==true){
        return values[4*t];
    }
    
    else{
        if(4*t+N_bins*4>8*N_bins-1){
            std::cout << "Warning: p0 exceeded the end of the density array, attempting to use index " << 4*t+N_bins*4 << ", t=" << t << std::endl;
        }
        return values[4*t+N_bins*4];
    }
}

void density::p_vector(int t, bool neutrino, three_vector* p)
{
    if(neutrino==true){
        for(int i=0; i<3; i++){
            p->set_value(i, values[4*t+i+1]);
        }}
    else{
        for(int i=0; i<3; i++){
            if(N_bins*4+4*t+1+i>8*N_bins-1){
                std::cout << "Warning: p_vector exceeded the end of the density array" << std::endl;
            }
            p->set_value(i, values[N_bins*4+4*t+1+i]);
        }}
}

void density::p0_p(int t, bool neutrino, three_vector* p)
{
    if(neutrino==true){
        for(int i=0; i<3; i++){
            p->set_value(i,values[4*t+i+1]);
        }
        p->multiply_by(values[4*t]);
        
    }
    
    else{
        for(int i=0; i<3; i++){
            p->set_value(i, values[N_bins*4+4*t+i+1]);
        }
        p->multiply_by(values[4*t+N_bins*4]);
    }
}

void density::number_density(double* output)
{
    dep_vars* nu_e = new dep_vars(N_bins);
    dep_vars* nu_mu = new dep_vars(N_bins);
    dep_vars* nubar_e = new dep_vars(N_bins);
    dep_vars* nubar_mu = new dep_vars(N_bins);

    double P0, P0bar, Pz, Pzbar, eps;
    for(int i = 0; i < N_bins; i++)
        {
            P0 = values[4*i];
            P0bar = values[4*i+N_bins*4];
            Pz = values[4*i+3];
            Pzbar = values[N_bins*4+4*i+3];
            eps = E->get_value(i);
            nu_e->set_value(i, 0.5 * P0 * (1 + Pz) * eps * eps);
            nu_mu->set_value(i, 0.5 * P0 * (1 - Pz) * eps * eps);
            nubar_e->set_value(i, 0.5 * P0bar * (1 + Pzbar) * eps * eps);
            nubar_mu->set_value(i, 0.5 * P0bar * (1 - Pzbar) * eps * eps);
        }

    double norm = pow(values[N-1], 3) / (2 * _PI_ * _PI_);
    output[0] = E->integrate(nu_e) * norm;
    output[1] = E->integrate(nu_mu) * norm;
    output[2] = E->integrate(nubar_e) * norm;
    output[3] = E->integrate(nubar_mu) * norm;

    delete nu_e;
    delete nu_mu;
    delete nubar_e;
    delete nubar_mu;
}

void density::print_csv(ostream& os)
{
    double nd[4];
    number_density(nd);

    for(int i = 0; i < 3; i++)
        os << nd[i] << ", ";
    os << nd[3];
}

integration::integration(linspace_and_gl* e, int p1_index){
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
            //count will give the number of energy values in eps that are less than or equal to the energy of p1+p2
            //this meants that count-1 will give the index of greatest element of eps less than the energy of p1+p2
            //furthermore, if count is bigger than N, the energy of p1+p2 is bigger than the biggest element of eps
            count = 0;
            for(int i=0; i<eps->get_len(); i++){
                if(eps->get_value(i)<=eps->get_value(p2)+eps->get_value(p1)){
                    count++;
                }
            }

            //p3_vals[p2] will have count+1 elements because we want count elements from eps as well as p1+p2 
            p3_vals[p2] = new dummy_vars(count+1);
            for(int i=0; i<count; i++){
                p3_vals[p2]->set_value(i, eps->get_value(i));
            }
            p3_vals[p2]->set_value(count,eps->get_value(p2)+eps->get_value(p1));
            p3_vals[p2]->set_trap_weights();

            inner_vals[p2] = new dep_vars(count+1);
        }
    }
    
    F_values = new double**[4]; 
    for(int i=0; i<4; i++){
        F_values[i] = new double*[eps->get_len()];
        for(int j=0; j<eps->get_len(); j++){
            F_values[i][j] = new double[eps->get_len()]();   
        }
    }
    
}


void integration::Fvvsc_components_term_1(density* dens, bool neutrino, int p2, int p3, double* F0, three_vector* F){
    
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();
    matrix* p_4 = new matrix();
    p_1->convert_p_to_identity_minus_matrix(dens, neutrino, p1);
    p_2->convert_p_to_identity_minus_matrix(dens, neutrino, p2);
    p_3->convert_p_to_matrix(dens, neutrino, p3);
    
    if(p1+p2-p3>=0){
        double p4_energy = eps->get_value(p1)+eps->get_value(p2)-eps->get_value(p3);
        //this clause finds an interpolated value for the p4 matrix if p4_energy is bigger than the biggest energy in the linspace
        double max_lin = eps->get_max_linspace();
        if (eps->get_value(p1)<=max_lin and eps->get_value(p2)<=max_lin and eps->get_value(p3)<=max_lin and p4_energy<=max_lin){
                p_4->convert_p_to_matrix(dens, neutrino, p1+p2-p3);
        }
        else{
            count = inner_vals[p2]->length()-1;
            p_4->convert_p4_to_interpolated_matrix(dens, neutrino, p4_energy, count);
        }
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


void integration::Fvvsc_components_term_2(density* dens, bool neutrino, int p2, int p3, double* F0, three_vector* F){
    
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();
    
    matrix* p_4 = new matrix(true);
    p_1->convert_p_to_matrix(dens, neutrino, p1);
    p_2->convert_p_to_matrix(dens, neutrino, p2);
    p_3->convert_p_to_identity_minus_matrix(dens, neutrino, p3);
    
    //this clause makes the p4 matrix the identity if p4 is bigger than max energy
    if(p1+p2-p3>=0){
        double p4_energy = eps->get_value(p1)+eps->get_value(p2)-eps->get_value(p3);
        //this clause finds an interpolated value for the p4 matrix if p4_energy is bigger than the biggest energy in the linspace
        double max_lin = eps->get_max_linspace();
        if (eps->get_value(p1)<=max_lin and eps->get_value(p2)<=max_lin and eps->get_value(p3)<=max_lin and p4_energy<=max_lin){
            p_4->convert_p_to_identity_minus_matrix(dens, neutrino, p1+p2-p3);
        }
        else{
            count = inner_vals[p2]->length()-1;
            p_4->convert_p4_to_identity_minus_interpolated_matrix(dens, neutrino, p4_energy, count);
        }
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

void integration::Fvvsc_components(density* dens, bool neutrino, int p2, int p3, double* F03, three_vector* F3){
    
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

void integration::all_F_for_p1(density* dens, bool neutrino){
    double F0 = 0;
    three_vector* Fxyz = new three_vector();
    for(int p2=0; p2<eps->get_len(); p2++){
        for(int p3=0; p3<eps->get_len(); p3++){
            
            if(p1+p2-p3>=0){
                Fvvsc_components(dens, neutrino, p2, p3, &F0, Fxyz);
                
                F_values[0][p2][p3] = F0;
                F_values[1][p2][p3] = Fxyz->get_value(0);
                F_values[2][p2][p3] = Fxyz->get_value(1);
                F_values[3][p2][p3] = Fxyz->get_value(2);
            }
        }
    }  
    for(int i=0; i<4; i++){
        for(int j=0; j<eps->get_len(); j++){
            for(int k=0; k<eps->get_len(); k++){
                if(isnan(F_values[i][j][k])){
                    std::cout << "F_values[" << i << "][" << j << "][" << k << "]=" << F_values[i][j][k] << std::endl;
                }}}}
    
    delete Fxyz;
}


double integration::J1(double p1, double p2, double p3){
    return 16./15 * pow(p3, 3) * (10 * pow(p1+p2, 2) - 15 * (p1+p2) * p3 + 6*pow(p3, 2));  
}

double integration::J2(double p1, double p2){
    return 16./15 * pow(p2, 3) * (10 * pow(p1, 2) + 5 * p1 * p2 + pow(p2, 2));  
}

double integration::J3(double p1, double p2, double p3){
    return 16./15 * (pow(p1+p2, 5) - 10 * pow(p1+p2, 2) * pow(p3, 3) + 15 * (p1+p2) * pow(p3, 4) - 6 * pow(p3,5));
}


double integration::interior_integral(density* dens, bool neutrino, int p2, int which_term){
    double p_1_energy = eps->get_value(p1);
    double max_energy = eps->get_value(eps->get_len()-1);
    

    //ALL OF THIS NEEDS TO BE TESTED AND DOUBLE CHECKED!!!!!!
    if(p2<p1){
        for(int p3=0; p3<p2; p3++){
            inner_vals[p2]->set_value(p3, F_values[which_term][p2][p3] * J1(p_1_energy, eps->get_value(p2), eps->get_value(p3)));
        }
        for(int p3=p2; p3<p1; p3++){
            inner_vals[p2]->set_value(p3, F_values[which_term][p2][p3] * J2(p_1_energy, eps->get_value(p2)));
        }
        
        //if p3 stays within the end of the linspace array
        if(eps->get_value(p2)+p_1_energy <= eps->get_max_linspace()){
            for(int p3=p1; p3<=p1+p2; p3++){
                inner_vals[p2]->set_value(p3, F_values[which_term][p2][p3] * J3(p_1_energy, eps->get_value(p2), eps->get_value(p3)));
            } 
        }
        else{
            for(int p3=p1; p3<inner_vals[p2]->length(); p3++){
                inner_vals[p2]->set_value(p3, F_values[which_term][p2][p3] * J3(p_1_energy, eps->get_value(p2), eps->get_value(p3)));
            }
        }
    }
    
    else{
        for(int p3=0; p3<p1; p3++){
            inner_vals[p2]->set_value(p3, F_values[which_term][p2][p3] * J1(p_1_energy, eps->get_value(p2), eps->get_value(p3)));
        }
        for(int p3=p1; p3<p2; p3++){
            inner_vals[p2]->set_value(p3, F_values[which_term][p2][p3] * J2(eps->get_value(p2), p_1_energy));
        }
        //if p3 stays within the end of the linspace array
        if(eps->get_value(p2)+p_1_energy <= eps->get_max_linspace()){
            for(int p3=p2; p3<=p1+p2; p3++){
                inner_vals[p2]->set_value(p3, F_values[which_term][p2][p3] * J3(p_1_energy, eps->get_value(p2), eps->get_value(p3)));
            } 
        }
        else{
            for(int p3=p2; p3<inner_vals[p2]->length(); p3++){
                inner_vals[p2]->set_value(p3, F_values[which_term][p2][p3] * J3(p_1_energy, eps->get_value(p2), eps->get_value(p3)));
            }
        }
    }
    double result = p3_vals[p2]->integrate(inner_vals[p2]);
    return result;
    
   
    
}

//note: results must be length 4
void integration::whole_integral(density* dens, bool neutrino, double* results){
    if (p1==0){
        for(int i=0; i<4; i++){
            results[i] = 0;
        }
    }
    else{
        //populates F_values
        all_F_for_p1(dens, neutrino);
        double Tcm = dens->get_Tcm();
            
        double p_1_energy = eps->get_value(p1);
        for(int i=0; i<4; i++){
            for(int p2=0; p2<eps->get_len(); p2++){
                outer_vals->set_value(p2, interior_integral(dens, neutrino, p2, i));
            }
            results[i] = eps->integrate(outer_vals);
            results[i] *= pow(Tcm, 5) * pow(_GF_,2) / (pow(2*_PI_,3) * pow(p_1_energy,2));

        }
       
    }
    
    
}


integration::~integration(){
    delete outer_vals;
    for(int i=0; i<eps->get_len(); i++){
        delete inner_vals[i];
        delete p3_vals[i];
    }
    delete[] inner_vals;
    delete[] p3_vals;
    for(int i=0; i<4; i++){
        for(int j=0; j<eps->get_len(); j++){
            delete[] F_values[i][j];
        }
        delete[] F_values[i];
    }
    delete[] F_values;
    delete eps;
    
}

/*
void density::print_csv(ostream& os)
{
    double nd[4];
    number_density(nd);

    three_vector_for_QKE* v = new three_vector_for_QKE;
    three_vector_for_QKE* v_vac = new three_vector_for_QKE;

    three_vector* p = new three_vector;
>>>>>>> chadbranch
    
    v->v_density(E, this);
    v_vac->v_vacuum(2.5e-15, 0.6, 0.8);
    v->add_to(-1., v_vac);

    p_vector(80, false, p);

    for(int i = 0; i < 3; i++)
        os << nd[i] << ", ";
    os << nd[3] << ", " << v->magnitude() << ", " << v->get_value(2) / v->magnitude() << ", ";
    p->print_csv(os);

    delete v;
    delete v_vac;
    delete p;
}
*/