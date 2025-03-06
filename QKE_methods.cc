#include "QKE_methods.hh"
#include <iostream>
#include "constants.hh"
#include <cmath>
#include "thermodynamics.hh"
#include "gl_vals.hh"
#include "matrices.hh"
#include <complex>

double extrapolate_exponential(double, double, double, double, double);
double extrapolate_linear(double, double, double, double, double);
double interpolate(double, int, double*, double*);


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
        values[4*i+3] =  (fnu - fmu)/(fnu+fmu+1.e-240);

        fnubar = 1 / (exp(eps_temp + eta_nu)+1);
        fmubar = 1 / (exp(eps_temp + eta_mu)+1);
        values[4*N_bins + 4*i] = fnubar + fmubar;
        values[4*N_bins + 4*i+3] = (fnubar - fmubar)/(fnu+fmu+1.e-240);
    }

}

density::density(dummy_vars* eps, int A, int B):dep_vars(8*eps->get_len()+2)
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
        fnu = (double)(A) * exp(-1 * eps_temp) / 10.;
        fmu = (double)(B) * exp(-2 * eps_temp) / 10.;
        values[4*i] = fnu + fmu;
        values[4*i+3] =  (fnu - fmu)/(fnu+fmu+1.e-240);

        fnubar = fnu;
        fmubar = fmu;
        values[4*N_bins + 4*i] = fnubar + fmubar;
        values[4*N_bins + 4*i+3] = (fnubar - fmubar)/(fnu+fmu+1.e-240);
    }

}

density::density(density* copy_me):dep_vars(copy_me)
{
    N_bins = copy_me->num_bins();
    E = new dummy_vars(copy_me->get_E());
}

density::~density(){    
    delete E; 
}

dummy_vars* density::get_E(){
    return E;
}

double density::get_E_value(int i){
    return E->get_value(i);
}

double density::get_T(){
    return values[N-2];
}

double density::get_Tcm(){
    return values[N-1];
}

int density::num_bins(){
    return N_bins;
}

void density::set_T(double T){ 
    values[N-2] = T;
    values[N-1] = T;
}

bool density::isnan(){
   for(int i=0; i< N; i++){
       if(std::isnan(values[i])){
           std::cout << "ERROR: density object is nan at index " << i << std::endl;
           return true;
       }
   }
    return false;
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

void density::p_vector(int t, bool neutrino, three_vector* p){
    if(neutrino==true){
        for(int i=0; i<3; i++){
            p->set_value(i, values[4*t+i+1]);
        }
    }
    else{
        for(int i=0; i<3; i++){
            if(N_bins*4+4*t+1+i>8*N_bins-1){
                std::cout << "Warning: p_vector exceeded the end of the density array" << std::endl;
            }
            p->set_value(i, values[N_bins*4+4*t+1+i]);
        }
    }
}

void density::p0_p(int t, bool neutrino, three_vector* p){
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

double density::interpolate_p0(bool neutrino, double energy){
    double interpolated_p0;

    int count = E->index_below_for_interpolation(energy);

    int back = 2;
    if (count - back < 0){
        back = count;
    }

    if(energy <= E->get_max_linspace()){
        double eps_values[4];
        double p_values[4];

        if (count == this->num_bins()-1){
            back = 3;
        }

        for(int j = 0; j < 4; j++){
            eps_values[j] = E->get_value(count-back+j);
            p_values[j] = std::log(this->p0(count-back+j, neutrino));
        }
        
        interpolated_p0 = std::exp(interpolate(energy, 4, eps_values, p_values));
    }
    
    else{
    //fixes indexing if count is last val in eps
        if(count==E->get_len()-1){
            count=E->get_len()-2;
        }
        interpolated_p0 = extrapolate_exponential(energy, E->get_value(count), E->get_value(count+1), this->p0(count, neutrino), this->p0(count+1, neutrino));
    }

    return interpolated_p0;
}

void density::interpolate_p0p(bool neutrino, double energy, three_vector* interpolated_p){
    three_vector** p_interp = new three_vector*[4];
    for(int j = 0; j < 4; j++)
    p_interp[j] = new three_vector();

    double temp_result;

    int count = E->index_below_for_interpolation(energy);


    int back = 2;
    if (count - back < 0){
        back = count;
    }
        
    //if p4 energy is below the max linspace
    if(energy <= E->get_max_linspace()){
        double p_values[4];
        double eps_values[4];

        if (count == this->num_bins()-1){
            back = 3;
        }


        for(int j=0; j<4; j++){
            eps_values[j] = E->get_value(count-back+j);
            p_values[j] = std::log(this->p0(count-back+j, neutrino));
            this->p0_p(count-back+j, neutrino, p_interp[j]);
            
            for(int i=0; i<3; i++){
                p_values[j] = p_interp[j]->get_value(i);
                temp_result = interpolate(energy, 4, eps_values, p_values);
                interpolated_p->set_value(i, temp_result);
            }
        }
    }
    else{
        //fixes indexing if count is last val in eps
        if(count==E->get_len()-1){
            count=E->get_len()-2;
        }

        this->p0_p(count, neutrino, p_interp[0]);
        this->p0_p(count+1, neutrino, p_interp[1]);
        for(int i=0; i<3; i++){
            temp_result = extrapolate_exponential(energy, E->get_value(count), E->get_value(count+1), p_interp[0]->get_value(i), p_interp[1]->get_value(i));
            interpolated_p->set_value(i, temp_result);
        }
    }

    for(int j = 0; j < 4; j++){
        delete p_interp[j];
    }
    delete[] p_interp;
}

double fifth_order_fit(double x, double* x_vals, double* y_vals){
    double fit = 0;
    double Lj;

    for(int j=0; j<5; j++){
        Lj = 1.0;
        for(int i=0; i<5; i++){
            if(i != j){
                Lj *= (x - x_vals[i]) / (x_vals[j] - x_vals[i]);
            }
        }
        fit += y_vals[j] * Lj;
    }
    return fit;
}

double interpolate_log_fifth(double x, double* x_vals, double* y_vals){
    double y_temp;


    for(int i=1; i<5; i++){
        if(y_vals[0] * y_vals[i] <= 0){
        return fifth_order_fit(x, x_vals, y_vals);
        }
    }

    double* y_log = new double[5]();
    for(int j=0; j<5; j++){
        y_log[j] = log(std::abs(y_vals[j]));
    }
    y_temp = fifth_order_fit(x, x_vals, y_log);
    
    delete[] y_log;

    if(y_vals[0] > 0){
        return exp(y_temp);
    }
    else{
        return -exp(y_temp);
    }

}


double linear(double x, double x1, double x2, double y1, double y2){
    if(x2-x1==0){std::cout << "warning: attempting to divide by 0**" << x << std::endl;}
        double slope = (y2-y1)/(x2-x1);
        return slope * (x-x2) + y2;
    }

double interpolate_log_linear(double x, double x_val1, double x_val2, double y_val1, double y_val2){

    if(y_val1 * y_val2 <= 0){
        return extrapolate_linear(x, x_val1, x_val2, y_val1, y_val2);
    }
    else{
        double y_temp = linear(x, x_val1, x_val2, log(std::abs(y_val1)), log(std::abs(y_val2)));
        if(y_val1 > 0){
            return exp(y_temp);
        }
        else{
            return -exp(y_temp);
        }
    }
}

double density::interpolated_matrix(bool neutrino, int index, double p4_energy, three_vector* p0p){
    double* results = new double[4]();
    double p0;

    //if in linspace, do fifth order interpolation
    if(p4_energy <= E->get_max_linspace()){
        int ind = std::max(0, index-2);
        ind = std::min(E->get_len()-1-4, ind);
        double* eps_vals = new double[5]();
        double* matrix_vals= new double[5]();

        for(int i=0; i<5; i++){
            eps_vals[i] = E->get_value(ind+i);
        }

        //1/2(p0+p0pz)
        for(int j=0; j<5; j++){
            this->p0_p(ind+j, neutrino, p0p);
            matrix_vals[j] = 0.5 * (this->p0(ind+j, neutrino) + p0p->get_value(2));
        }
        results[0] = interpolate_log_fifth(p4_energy, eps_vals, matrix_vals);

        //1/2p0px
        for(int j=0; j<5; j++){
            this->p0_p(ind+j, neutrino, p0p);
            matrix_vals[j] = 0.5 * p0p->get_value(0);
        }
        results[1] = interpolate_log_fifth(p4_energy, eps_vals, matrix_vals);

        //1/2p0py
        for(int j=0; j<5; j++){
            this->p0_p(ind+j, neutrino, p0p);
            matrix_vals[j] = 0.5 * p0p->get_value(1);
        }
        results[2] = interpolate_log_fifth(p4_energy, eps_vals, matrix_vals);

        //1/2(p0-p0pz)
        for(int j=0; j<5; j++){
            this->p0_p(ind+j, neutrino, p0p);
            matrix_vals[j] = 0.5 * (this->p0(ind+j, neutrino) - p0p->get_value(2));
        }
        results[3] = interpolate_log_fifth(p4_energy, eps_vals, matrix_vals);


        delete[] eps_vals;
        delete[] matrix_vals;
    }


    //if not in linspace do linear interpolation
    else{
        if(index==E->get_len()-1){
            index = index-1;
        }
        double energy_one = E->get_value(index);
        double energy_two = E->get_value(index+1);
        three_vector* secondp0p = new three_vector();

        this->p0_p(index, neutrino, p0p);
        this->p0_p(index+1, neutrino, secondp0p);

        //1/2(p0+p0pz)
        results[0] = interpolate_log_linear(p4_energy, energy_one, energy_two, 0.5*(this->p0(index, neutrino) + p0p->get_value(2)), 0.5*(this->p0(index+1, neutrino) + secondp0p->get_value(2)));

        //1/2(p0px)
        results[1] = interpolate_log_linear(p4_energy, energy_one, energy_two, 0.5*p0p->get_value(0), 0.5*secondp0p->get_value(0));

        //1/2(p0py)
        results[2] = interpolate_log_linear(p4_energy, energy_one, energy_two, 0.5*p0p->get_value(1), 0.5*secondp0p->get_value(1));

        //1/2(p0-p0pz)
        results[3] = interpolate_log_linear(p4_energy, energy_one, energy_two, 0.5*(this->p0(index, neutrino) - p0p->get_value(2)), 0.5*(this->p0(index+1, neutrino) - secondp0p->get_value(2)));

        delete secondp0p;
    }


    p0p->set_value(0, 2*results[1]);
    p0p->set_value(1, 2*results[2]);
    p0p->set_value(2, results[0]-results[3]);

    
    //p0
    p0 = results[0] + results[3];
    delete[] results;

    
    return p0;

}

double interpolate(double x, int N, double* x_vals, double* y_vals){
    double res = 0;
    double termj = 1;
    for (int j = 0; j < N; j++)
    {
        termj = 1;
        for(int k = 0; k < N; k++){
            if (j!=k){
                termj *= (x - x_vals[k]) / (x_vals[j] - x_vals[k]);
                termj *= y_vals[j];
                res += termj;
            }
        }
    }
    return res;
}

double extrapolate_exponential(double x, double x1, double x2, double y1, double y2){
    //note: this assumes x1<x2, so we expect y1>y2 because this is an exponential decay model
    if(y1==y2){
        return y1;
    }

    else{
        //model is Ce^(-ax)
        if(y1/y2 < 1){
            return extrapolate_linear(x, x1, x2, y1, y2);
        }
        else{
            if(x1-x2==0){std::cout << "warning: attempting to divide by 0" << x << std::endl;}
            double a = -log(y1/y2) / (x1-x2);
            double C = y1 * exp(a * x1);
            return C * exp(-a * x);
        }
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



nu_nu_collision::nu_nu_collision(linspace_and_gl* e, int p1_index){
    eps = new linspace_and_gl(e);
    //eps = new linspace_and_gl_booles(e);
    p1 = p1_index;
    double p1_energy = eps->get_value(p1);
    outer_vals = new dep_vars(eps->get_len());
    inner_vals = new dep_vars*[eps->get_len()];
    p3_vals = new dummy_vars*[eps->get_len()];
    
    for(int p2=0; p2<eps->get_len(); p2++){
        
        p3_vals[p2] = new dummy_vars(eps);
        inner_vals[p2] = new dep_vars(eps->get_len());
        /*
        p3_vals[p2] = new gel_dummy_vars(100, eps->get_value(0), p1_energy+eps->get_value(p2));
        inner_vals[p2] = new dep_vars(p3_vals[p2]->get_len());
        */
    }
    
    //interpolation_indices[p2][p3]=[p3 index for interpolation, p4 index for interpolation]
    interpolation_indices = new int**[eps->get_len()];
    for(int p2=0; p2<eps->get_len(); p2++){
        interpolation_indices[p2] = new int*[p3_vals[p2]->get_len()];
        for(int p3=0; p3<p3_vals[p2]->get_len(); p3++){
            interpolation_indices[p2][p3] = new int[2];
        }
    }
    
    double p2_energy;
    double p3_energy;
    double p4_energy;
    
    for(int p2=0; p2<eps->get_len(); p2++){
        p2_energy = eps->get_value(p2);
        for(int p3=0; p3<p3_vals[p2]->get_len(); p3++){
            p3_energy = p3_vals[p2]->get_value(p3);
            p4_energy = p1_energy + p2_energy - p3_energy;
            interpolation_indices[p2][p3][0] = eps->index_below_for_interpolation(p3_energy);
            if(p4_energy>=0){
                interpolation_indices[p2][p3][1] = eps->index_below_for_interpolation(p4_energy);            
            }
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


void nu_nu_collision::Fvvsc_components_term_1(density* dens, bool neutrino, int p2, int p3, double* F0, three_vector* F){
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();
    matrix* p_4 = new matrix();
    p_1->convert_p_to_identity_minus_matrix(dens, neutrino, p1);
    p_2->convert_p_to_identity_minus_matrix(dens, neutrino, p2);

    double p3_energy = p3_vals[p2]->get_value(p3);
    double p4_energy = eps->get_value(p1)+eps->get_value(p2)-p3_energy;
    if(p4_energy<0){
        p4_energy = 0;
    }
    
    three_vector* A = new three_vector();
    double A0 = dens->interpolated_matrix(neutrino, interpolation_indices[p2][p3][0], p3_energy, A);
    p_3->convert_p_to_matrix(A0,A);

    A0 = dens->interpolated_matrix(neutrino, interpolation_indices[p2][p3][1], p4_energy, A);
    p_4->convert_p_to_matrix(A0,A);
    delete A;

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


void nu_nu_collision::Fvvsc_components_term_2(density* dens, bool neutrino, int p2, int p3, double* F0, three_vector* F){
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();

    matrix* p_4 = new matrix(true);
    p_1->convert_p_to_matrix(dens, neutrino, p1);
    p_2->convert_p_to_matrix(dens, neutrino, p2);
    
    double p3_energy = p3_vals[p2]->get_value(p3);
    double p4_energy = eps->get_value(p1)+eps->get_value(p2)-p3_energy;
    if(p4_energy<0){
        p4_energy = 0;
    }
    
    three_vector* A = new three_vector();
    double A0 = dens->interpolated_matrix(neutrino, interpolation_indices[p2][p3][0], p3_energy, A);
    p_3->convert_p_to_identity_minus_matrix(A0,A);

    A0 = dens->interpolated_matrix(neutrino, interpolation_indices[p2][p3][1], p4_energy, A);
    p_4->convert_p_to_identity_minus_matrix(A0,A);
    delete A;

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

void nu_nu_collision::Fvvsc_components(density* dens, bool neutrino, int p2, int p3, double* F03, three_vector* F3, bool net){
    double F01;
    three_vector* F1 = new three_vector();
    double F02;
    three_vector* F2 = new three_vector();


    Fvvsc_components_term_1(dens, neutrino, p2, p3, &F01, F1);
    Fvvsc_components_term_2(dens, neutrino, p2, p3, &F02, F2);

    if(net==true){
        F2->multiply_by(-1);
        *F03 = F01 - F02;
    }
    else{
        *F03 = F01 + F02;
    }
    F3->add(F1, F2);

    delete F1;
    delete F2;
}

void nu_nu_collision::Fvvsc_for_p1(density* dens, bool neutrino, bool net){
    double F0 = 0;
    three_vector* Fxyz = new three_vector();
    for(int p2=0; p2<eps->get_len(); p2++){
        for(int p3=0; p3<p3_vals[p2]->get_len(); p3++){
            //this clause means F is filled in only if p3_energy is less than p1_energy+p2_energy
            //this won't affect results for p3 objects that go up to p1+p2, and this is useful in trap rule when we want to set integrand to 0 past p1+p2
            if(p3_vals[p2]->get_value(p3) <= eps->get_value(p1) + eps->get_value(p2)){
                //if p4_energy>=0
                if(eps->get_value(p1)+eps->get_value(p2)-p3_vals[p2]->get_value(p3)>=0){
                    Fvvsc_components(dens, neutrino, p2, p3, &F0, Fxyz, net);

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

void nu_nu_collision::Fvvbarsc_components_term_1(density* dens, bool neutrino, int p2, int p3, double* F0, three_vector* F){
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();

    matrix* p_4 = new matrix(true);
    p_1->convert_p_to_identity_minus_matrix(dens, neutrino, p1);
    p_2->convert_p_to_identity_minus_matrix(dens, not neutrino, p2);
    
    double p3_energy = p3_vals[p2]->get_value(p3);
    double p4_energy = eps->get_value(p1)+eps->get_value(p2)-p3_energy;
    if(p4_energy<0){
        p4_energy = 0;
    }
    
    three_vector* A = new three_vector();
    double A0 = dens->interpolated_matrix(not neutrino, interpolation_indices[p2][p3][0], p3_energy, A);
    p_3->convert_p_to_matrix(A0,A);

    A0 = dens->interpolated_matrix(neutrino, interpolation_indices[p2][p3][1], p4_energy, A);
    p_4->convert_p_to_matrix(A0,A);
    delete A;

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

void nu_nu_collision::Fvvbarsc_components_term_2(density* dens, bool neutrino, int p2, int p3, double* F0, three_vector* F){
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();

    matrix* p_4 = new matrix(true);
    p_1->convert_p_to_matrix(dens, neutrino, p1);
    p_2->convert_p_to_matrix(dens, not neutrino, p2);
    
    double p3_energy = p3_vals[p2]->get_value(p3);
    double p4_energy = eps->get_value(p1)+eps->get_value(p2)-p3_energy;
    if(p4_energy<0){
        p4_energy = 0;
    }
    
    three_vector* A = new three_vector();
    double A0 = dens->interpolated_matrix(not neutrino, interpolation_indices[p2][p3][0], p3_energy, A);
    p_3->convert_p_to_identity_minus_matrix(A0,A);

    A0 = dens->interpolated_matrix(neutrino, interpolation_indices[p2][p3][1], p4_energy, A);
    p_4->convert_p_to_identity_minus_matrix(A0,A);
    delete A;

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

void nu_nu_collision::Fvvbarsc_components(density* dens, bool neutrino, int p2, int p3, double* F03, three_vector* F3, bool net){
    double F01;
    three_vector* F1 = new three_vector();
    double F02;
    three_vector* F2 = new three_vector();

    Fvvbarsc_components_term_1(dens, neutrino, p2, p3, &F01, F1);
    Fvvbarsc_components_term_2(dens, neutrino, p2, p3, &F02, F2);

    if(net==true){
        F2->multiply_by(-1);
        *F03 = F01 - F02;
    }
    else{
        *F03 = F01 + F02;
    }

    F3->add(F1, F2);


    delete F1;
    delete F2;
}

void nu_nu_collision::Fvvbarsc_for_p1(density* dens, bool neutrino, bool net){
    double F0 = 0;
    three_vector* Fxyz = new three_vector();
    for(int p2=0; p2<eps->get_len(); p2++){
        for(int p3=0; p3<p3_vals[p2]->get_len(); p3++){
            //only if p3_energy is less than p1_energy+p2_energy is this called--> will make integrand 0 past p1+p2 for trap 
            //won't affect results for other dummy var objects
            if(p3_vals[p2]->get_value(p3) <= eps->get_value(p1) + eps->get_value(p2)){
                //p4_energy must be >=0
                if(eps->get_value(p1)+eps->get_value(p2)-p3_vals[p2]->get_value(p3)>=0){
                    Fvvbarsc_components(dens, neutrino, p2, p3, &F0, Fxyz, net);

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


double nu_nu_collision::J1(double p1, double p2, double p3){
    return 16./15 * pow(p3,3) * (10 * pow(p1+p2,2) - 15 * (p1+p2) * p3 + 6*pow(p3,2));  
}

double nu_nu_collision::J2(double p1, double p2){
    return 16./15 * pow(p2,3) * (10 * pow(p1,2) + 5 * p1*p2 + pow(p2,2));  
}

double nu_nu_collision::J3(double p1, double p2, double p3){
    return 16./15 * (pow(p1+p2,5) - 10 * pow(p1+p2,2) * pow(p3, 3) + 15 * (p1+p2) * pow(p3,4) - 6 * pow(p3,5));
}

double nu_nu_collision::K1(double p1, double p3){
    return 16./15 * pow(p3,3) * (10 * pow(p1,2) - 5 * p1*p3 + pow(p3,2));
}

double nu_nu_collision::K2(double p1, double p2, double p3){
    return 16./15 * pow(p2,3) * (10 * pow(p1-p3,2) + 15 * (p1-p3) * p2 + 6 * pow(p2,2));
}

double nu_nu_collision::K3(double p1, double p2, double p3){
    return 16./15 * (pow(p1-p3,5) + 10 * pow(p1-p3,2) * pow(p2,3) + 15 * (p1-p3) * pow(p2,4) + 6 * pow(p2,5));
}


double nu_nu_collision::interior_integral(int p2, int which_term){
    double p_1_energy = eps->get_value(p1);
    double max_energy = eps->get_value(eps->get_len()-1);
    double p2_energy = eps->get_value(p2);
    
    int p3 = 0;
    double p3_energy;
    
    while(p3 < p3_vals[p2]->get_len()){
        p3_energy = p3_vals[p2]->get_value(p3);
    
        if(p2<p1){
            if(p3_energy < p2_energy){
                inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J1(p_1_energy, p2_energy, p3_energy) + Fvvbar_values[which_term][p2][p3] * K1(p_1_energy, p3_energy));
               
            }
            else if(p3_energy < eps->get_value(p1)){
                inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J2(p_1_energy, p2_energy) + Fvvbar_values[which_term][p2][p3] * K2(p_1_energy, p2_energy, p3_energy));
            }
            else{
                inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J3(p_1_energy, p2_energy, p3_energy) + Fvvbar_values[which_term][p2][p3] * K3(p_1_energy, p2_energy, p3_energy));
            }

        }

        else{
            if(p3_energy < p_1_energy){
                inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J1(p_1_energy, p2_energy, p3_energy) + Fvvbar_values[which_term][p2][p3] * K1(p_1_energy, p3_energy));

            }
            else if(p3_energy < p2_energy){
                inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J2(p2_energy, p_1_energy) + Fvvbar_values[which_term][p2][p3] * K1(p3_energy, p_1_energy));

            }
            else{
                inner_vals[p2]->set_value(p3, Fvv_values[which_term][p2][p3] * J3(p_1_energy, p2_energy, p3_energy) + Fvvbar_values[which_term][p2][p3] * K3(p_1_energy, p2_energy, p3_energy));
            }
        }
        p3++;
    }
        
    double result = p3_vals[p2]->integrate(inner_vals[p2]);
    return result;
}

//note: results must be length 4
void nu_nu_collision::whole_integral(density* dens, bool neutrino, double* results, bool net){
    if (p1==0){
        for(int i=0; i<4; i++){
            results[i] = 0;
        }
    }
    else{
        //populates F_values
        Fvvsc_for_p1(dens, neutrino, net);
        Fvvbarsc_for_p1(dens, neutrino, net);
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


nu_nu_collision::~nu_nu_collision(){
    for(int i=0; i<eps->get_len(); i++){
        for(int j=0; j<p3_vals[i]->get_len(); j++){
            delete[] interpolation_indices[i][j];
        }
        delete[] interpolation_indices[i];
    }
    delete[] interpolation_indices;
    
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



nu_e_collision::nu_e_collision(linspace_and_gl* e, int p1_index, double T_comoving){
   Tcm = T_comoving;
   scaled_me = _electron_mass_ / Tcm;
   me_squared = pow(scaled_me,2);
   eps = new linspace_and_gl(e);
   p1 = p1_index;
   p1_energy = eps->get_value(p1);
   p1_me = p1_energy / scaled_me;
   
   int numgl_points = 50;
   outer_vals_R2 = new dep_vars(numgl_points);
   q3_vals_R2 = new dummy_vars(numgl_points);
   outer_vals_R1 = new dep_vars(numgl_points);
   q2_vals_R1 = new dummy_vars(numgl_points);
   for(int i=0; i<numgl_points; i++){
       q3_vals_R2->set_value(i, xvals_50[i]);
       q3_vals_R2->set_weight(i, wvals_50[i]*exp(xvals_50[i]));
       q2_vals_R1->set_value(i, xvals_50[i]);
       q2_vals_R1->set_weight(i, wvals_50[i]*exp(xvals_50[i]));
   }
   
   double E_cut_1_R2 = p1_energy + me_squared/(4*p1_energy);
   double E_cut_2_R2 = p1_energy + scaled_me*(p1_energy+scaled_me)/(2*p1_energy+scaled_me);
   double E_cut_3 = sqrt(pow(p1_energy,2) + me_squared);
   q_cut_1_R2 = sqrt(pow(E_cut_1_R2,2) - me_squared);
   q_cut_2_R2 = sqrt(pow(E_cut_2_R2,2) - me_squared);
   q_cut_3 = sqrt(pow(E_cut_3,2) - me_squared);
   
   q_trans_2_R2 = new dep_vars(q3_vals_R2->get_len());
   q_lim_1_R2 = new dep_vars(q3_vals_R2->get_len());
   double q3_momentum = 0;
   double E3 = 0;
   double E_lim_1_R2 = 0;
   double E_trans_2_R2;
   
   for(int q3=0; q3<q3_vals_R2->get_len(); q3++){
       q3_momentum = q3_vals_R2->get_value(q3);
       E3 = sqrt(pow(q3_momentum,2) + me_squared);
       E_trans_2_R2 = 0.5 * (E3 + q3_momentum - 2*p1_energy + me_squared / (E3 + q3_momentum - 2*p1_energy));
       E_lim_1_R2 = 0.5 * (E3 - q3_momentum - 2*p1_energy + me_squared / (E3 - q3_momentum - 2*p1_energy));
       
       q_trans_2_R2->set_value(q3, sqrt(pow(E_trans_2_R2,2) - me_squared));
       q_lim_1_R2->set_value(q3, sqrt(pow(E_lim_1_R2,2) - me_squared));
   }
   
   double E_cut_1_R1 = scaled_me + 2*pow(p1_energy,2)/(scaled_me - 2*p1_energy);
   q_cut_1_R1 = sqrt(pow(E_cut_1_R1,2) + me_squared);
   q_trans_2_R1 = new dep_vars(q2_vals_R1->get_len());
   q_lim_1_R1 = new dep_vars(q2_vals_R1->get_len());
   double q2_momentum = 0;
   double E2 = 0;
   double E_lim_1_R1 = 0;
   double E_trans_2_R1;
   
   for(int q2=0; q2<q2_vals_R1->get_len(); q2++){
       q2_momentum = q2_vals_R1->get_value(q2);
       E2 = sqrt(pow(q2_momentum,2) + me_squared);
       E_trans_2_R1 = 0.5 * (2*p1_energy + E2 - q2_momentum + me_squared / (2*p1_energy + E2 - q2_momentum));
       E_lim_1_R1 = 0.5 * (2*p1_energy + E2 + q2_momentum + me_squared / (2*p1_energy + E2 + q2_momentum));
       
       q_trans_2_R1->set_value(q2, sqrt(pow(E_trans_2_R1,2) - me_squared));
       q_lim_1_R1->set_value(q2, sqrt(pow(E_lim_1_R1,2) - me_squared));
   }
   
   p4_min_vals_R2 = new double[q3_vals_R2->get_len()]();
   p4_max_vals_R2 = new double[q3_vals_R2->get_len()]();
   count_min_vals_R2 = new int[q3_vals_R2->get_len()]();
   count_max_vals_R2 = new int[q3_vals_R2->get_len()]();
   
   double p4_min = 0;
   double p4_max = 0;
   int count_min = 0;
   int count_max = 0;
   double q3_energy = 0;
   double q2_min = 0;
   double q2_max = 0;
   
   inner_vals_R2 = new dep_vars*[outer_vals_R2->length()];
   q2_vals_R2 = new dummy_vars*[q3_vals_R2->get_len()];
   
   
   for(int q3=0; q3<q3_vals_R2->get_len(); q3++){
       q3_momentum = q3_vals_R2->get_value(q3);
       E3 = sqrt(pow(q3_momentum,2) + me_squared);
       //idea here is to establish what we want p4 vals to be and then use those vals to reconstruct E3 vals
       
       //first we decide minimum and maximum p4 values. these will be only interpolated p4 values
       //we have to consider cases 
 
       //case 1
       if(p1_me < (sqrt(5)-1)/4){
           //case 1a: m_e < E3 < E_cut_2 => m_e < E2 < E_lim_1
           if(E3 < E_cut_2_R2){
               q2_min = 0;
               q2_max = q_lim_1_R2->get_value(q3);
           }
           else{
               q2_min = q_trans_2_R2->get_value(q3);
               //case 1b: E_cut_2 < E3 < E_cut_1 => E_lim_2 < E2 < E_lim_1                  
               if(E3 < E_cut_1_R2){
                   q2_max = q_lim_1_R2->get_value(q3);
               }
               //case 1c: E_cut_1 < E3 < inf => E_lim_2 < E2 < inf
               else{
                   q2_max = q3_vals_R2->get_value(numgl_points-1);
               } 
           }
       }
       //case 2
       else if(p1_me < 1/(2*sqrt(2))){
           if(E3 < E_cut_2_R2){
               q2_min = 0;
               
               //case 2a: m_e < E3 < E_cut_1 => m_e < E2 < E_lim_1
               if(E3 < E_cut_1_R2){
                   q2_max = q_lim_1_R2->get_value(q3);
               }
               //case 2b: m_e < E3 < inf => m_e < E2 < inf
               else{
                   q2_max = q3_vals_R2->get_value(numgl_points-1);
               }
               
           }
           //case 2c: E_cut_2 < E3 < inf => E_lim_2 < E2 < inf
           else{
               q2_min = q_trans_2_R2->get_value(q3);
               q2_max = q3_vals_R2->get_value(numgl_points-1);
           }
       }
       
       //case 3
       else if(p1_me < 0.5){
           if(E3 < E_cut_2_R2){   
               q2_min = 0;
               //case 3a: m_e < E3 < E_cut_1 => m_e < E2 < E_lim_1
               if(E3 < E_cut_1_R2){
                   q2_max = q_lim_1_R2->get_value(q3);
               }
               //case 3b: E_cut_1 < E3 < E_cut_2 => m_e < E2 < inf
               else{
                   q2_max = q3_vals_R2->get_value(numgl_points-1);
               }
           }
           //case 3c: E_cut_2 < E3 < inf => E_lim_2 < E2 < inf
           else{
               q2_min = q_trans_2_R2->get_value(q3);
               q2_max = q3_vals_R2->get_value(numgl_points-1);
           }
       }
       //case 4
       else{
           q2_max = q3_vals_R2->get_value(numgl_points-1);
           //case 4a: m_e < E3 < E_cut_2 => m_e < E2 < inf
           if(E3 < E_cut_2_R2){
               q2_min = 0;
           }
           //case 4b: E_cut_2 < E3 < inf => E_lim_2 < E2 < inf
           else{
               q2_min = q_trans_2_R2->get_value(q3);
           }
       }
       
       p4_min = p1_energy + sqrt(pow(q2_min,2) + me_squared) - sqrt(pow(q3_momentum,2) + me_squared);
       p4_max = p1_energy + sqrt(pow(q2_max,2) + me_squared) - sqrt(pow(q3_momentum,2) + me_squared);
       
       //WARNING: potential issues in the code that result in p4_min or p4_max being negative WILL NOT BE CAUGHT
       if(p4_min < 0){
           if(p4_min < -1){
               std::cout << "warning: p4_min is " << p4_min << ", setting to 0" << std::endl;
           }
           p4_min=0;
       }
       if(p4_max < 0){
           if(p4_max < -1){
               std::cout << "warning: p4_max is " << p4_max << ", setting to 0" << std::endl;
           }
           p4_max=0;
       }
       
       p4_min_vals_R2[q3] = p4_min;
       p4_max_vals_R2[q3] = p4_max;
       
       double temp_energy = eps->get_value(0);
       count_min = eps->index_below_for_interpolation(p4_min)+1;
       count_max = eps->index_below_for_interpolation(p4_max);
       //count_min gives the number of items in epsilon that have energy less than the minimum p4 val; therefore first p4 val of interest is epsilon[count_min]
       //count_max gives the index of the greatest element of epsilon that has energy less than the maximum p4 val; therefore last p4 val of interest is epsilon[count_max]
       
       //p4 vals will contain p4_min, epsilon values from indices count_min to count_max, inclusive, and p4_max
       //therefore E2_vals needs to have 3+count_max-count_min things in it
       //note that if count_min and count_max beyond the end of the array, count_max=count_min-1 so q2_vals_R2 will have length 2
       
       count_min_vals_R2[q3] = count_min;
       count_max_vals_R2[q3] = count_max;
       
       q2_vals_R2[q3] = new dummy_vars(count_max-count_min+3);
       q2_vals_R2[q3]->set_value(0, q2_min);
       q2_vals_R2[q3]->set_value(count_max-count_min+2, q2_max);
       
       for(int j=count_min; j<=count_max; j++){
           //q2 = p4 - p1 + q3
           q2_vals_R2[q3]->set_value(j-count_min+1, sqrt(pow(eps->get_value(j) - p1_energy + E3,2) - me_squared));
       }
       
       q2_vals_R2[q3]->set_trap_weights();
       inner_vals_R2[q3] = new dep_vars(count_max-count_min+3);
   }
   
   inner_vals_R1 = new dep_vars*[outer_vals_R1->length()];
   q3_vals_R1 = new dummy_vars*[q2_vals_R1->get_len()];
   
   p4_min_vals_R1 = new double[q2_vals_R1->get_len()]();
   p4_max_vals_R1 = new double[q2_vals_R1->get_len()]();
   count_min_vals_R1 = new int[q2_vals_R1->get_len()]();
   count_max_vals_R1 = new int[q2_vals_R1->get_len()]();
   
   double q2_energy = 0;
   double q3_min = 0;
   double q3_max = 0;
   
   for(int q2=0; q2<q2_vals_R1->get_len(); q2++){
       q2_momentum = q2_vals_R1->get_value(q2);
       E2 = sqrt(pow(q2_momentum,2) + me_squared);
       
       q3_max = q_lim_1_R1->get_value(q2);
       if(p1_me < 0.5 and E2 > E_cut_1_R1){
           //case 1b: E_cut_1 < E2 < inf => E_lim_1 < E3 < E_lim_1
           q3_min = 0;
       } 
       else{
           //case 1a: m_e < E2 < E_cut_1 => m_e < E3 < E_lim_1
           //case 2
           q3_min = q_trans_2_R1->get_value(q2);
       }
       
       
       p4_min = p1_energy + E2 - sqrt(pow(q3_max,2) + me_squared);
       p4_max = p1_energy + E2 - sqrt(pow(q3_min,2) + me_squared);
       
       //WARNING: potential issues in the code that result in p4_min or p4_max being negative WILL NOT BE CAUGHT
       if(p4_min < 0){
           if(p4_min < -1){
               std::cout << "warning: p4_min is " << p4_min << ", setting to 0" << std::endl;
           }
           p4_min=0;
       }
       if(p4_max < 0){
           if(p4_max < -1){
               std::cout << "warning: p4_max is " << p4_max << ", setting to 0" << std::endl;
           }
           p4_max=0;
       }
       
       p4_min_vals_R1[q2] = p4_min;
       p4_max_vals_R1[q2] = p4_max;
       
       double temp_energy = eps->get_value(0);
       count_min = eps->index_below_for_interpolation(p4_min)+1;
       count_max = eps->index_below_for_interpolation(p4_max);
       //count_min gives the number of items in epsilon that have energy less than the minimum p4 val; therefore first p4 val of interest is epsilon[count_min]
       //count_max gives the index of the greatest element of epsilon that has energy less than the maximum p4 val; therefore last p4 val of interest is epsilon[count_max]
       
       //p4 vals will contain p4_min, epsilon values from indices count_min to count_max, inclusive, and p4_max
       //therefore E2_vals needs to have 3+count_max-count_min things in it
       //note that if count_min and count_max beyond the end of the array, count_max=count_min-1 so q2_vals_R2 will have length 2
       count_min_vals_R1[q2] = count_min;
       count_max_vals_R1[q2] = count_max;
       
       
       q3_vals_R1[q2] = new dummy_vars(count_max-count_min+3);
       q3_vals_R1[q2]->set_value(0, q3_min);
       q3_vals_R1[q2]->set_value(count_max-count_min+2, q3_max);
       
       for(int j=count_min; j<=count_max; j++){
           //q3 = p1 + q2 - p4
           q3_vals_R1[q2]->set_value(count_max-j+1, sqrt(pow(p1_energy + E2 - eps->get_value(j),2) - me_squared));
       }
       
       q3_vals_R1[q2]->set_trap_weights();
       inner_vals_R1[q2] = new dep_vars(count_max-count_min+3);
       
       
   }
   
   R2_F_LL_RR_values = new double**[4]();
   R2_F_LR_RL_values = new double**[4]();
   for(int i=0; i<4; i++){
       R2_F_LL_RR_values[i] = new double*[eps->get_len()+1]();
       R2_F_LR_RL_values[i] = new double*[eps->get_len()+1]();
       for(int j=0; j<eps->get_len()+1; j++){
           R2_F_LL_RR_values[i][j] = new double[q3_vals_R2->get_len()];
           R2_F_LR_RL_values[i][j] = new double[q3_vals_R2->get_len()];
       }
       
   }
                                     
                                     
   R1_F_LL_RR_values = new double**[4]();
   R1_F_LR_RL_values = new double**[4]();
   for(int i=0; i<4; i++){
       R1_F_LL_RR_values[i] = new double*[q2_vals_R1->get_len()]();
       R1_F_LR_RL_values[i] = new double*[q2_vals_R1->get_len()]();
       for(int j=0; j<q2_vals_R1->get_len(); j++){
           R1_F_LL_RR_values[i][j] = new double[eps->get_len()+1];
           R1_F_LR_RL_values[i][j] = new double[eps->get_len()+1];
       }
       
   }
   
}


void nu_e_collision::F_LL_F_RR(double* F0, three_vector* F, density* dens, bool neutrino, int q2, double E2, int q3, double E3, int p4, double p4_energy, int count_min, int count_max){
   complex_three_vector* A = new complex_three_vector();
   A->set_value(2, complex<double> (0.5,0));
   matrix* G_L = new matrix(complex<double> (_sin_squared_theta_W_,0),A);
   
   matrix* G_R = new matrix(true);
   G_R->multiply_by(_sin_squared_theta_W_);
   
   matrix* p_1 = new matrix();
   matrix* minus_p_1 = new matrix();
   matrix* p_4 = new matrix();
   matrix* minus_p_4 = new matrix(true);
   
   p_1->convert_p_to_matrix(dens, neutrino, p1);
   minus_p_1->convert_p_to_identity_minus_matrix(dens, neutrino, p1);
   
   double f2 = 1 / (exp(E2/Tcm)+1);
   double f3 = 1 / (exp(E3/Tcm)+1);
   /*
   //case of p4min
   if(p4 == -1){
       p_4->convert_p4_to_interpolated_matrix(dens, neutrino, p4_energy);
       minus_p_4->convert_p4_to_identity_minus_interpolated_matrix(dens, neutrino, p4_energy);
   }
   //case of p4max
   else if(p4 == -2){
       p_4->convert_p4_to_interpolated_matrix(dens, neutrino, p4_energy);
       minus_p_4->convert_p4_to_identity_minus_interpolated_matrix(dens, neutrino, p4_energy);
   }
   else{
       p_4->convert_p_to_matrix(dens, neutrino, p4);
       minus_p_4->convert_p_to_identity_minus_matrix(dens, neutrino, p4);
   }
   */
/*
   F_dummy1 = G_L * rho_4
   F_dummy2 = G_L * (1-rho_1)
   F_dummy3 = F_dummy1 * F_dummy2
   F_dummy4 = G_L * (1-rho_4)
   F_dummy5 = G_L * rho_1
   F_dummy6 = F_dummy4 * F_dummy5
   F_dummy7 = F_dummy3 - F_dummy6 => F_dummy7 = F_LL
   
   F_dummy8 = rho_4 * (1-rho_1)
   F_dummy9 = (1-rho_4) * rho_1
   F_dummy10 = F_dummy8 - F_dummy10 => F_dummy10 = F_RR
   
   F_dummy11 = F_dummy7 + F_dummy10 = F_LL + F_RR
   */

   
   matrix* F_dummy1 = new matrix();
   matrix* F_dummy2 = new matrix();
   matrix* F_dummy3 = new matrix();
   matrix* F_dummy4 = new matrix();
   matrix* F_dummy5 = new matrix();
   matrix* F_dummy6 = new matrix();
   matrix* F_dummy7 = new matrix();
   matrix* F_dummy8 = new matrix();
   matrix* F_dummy9 = new matrix();
   matrix* F_dummy10 = new matrix();
   matrix* F_dummy11 = new matrix();
   
   F_dummy1->matrix_multiply(G_L, p_4);
   F_dummy2->matrix_multiply(G_L, minus_p_1);
   F_dummy3->matrix_multiply(F_dummy1, F_dummy2);
   F_dummy3->multiply_by(f3 * (1-f2));
   F_dummy4->matrix_multiply(G_L, minus_p_4);
   F_dummy5->matrix_multiply(G_L, p_1);
   F_dummy6->matrix_multiply(F_dummy4, F_dummy5);
   F_dummy6->multiply_by(f2 * (1-f3));
   F_dummy6->multiply_by(complex<double> (-1,0));
   F_dummy7->matrix_add(F_dummy3, F_dummy6);
   
   F_dummy8->matrix_multiply(p_4, minus_p_1);
   F_dummy8->multiply_by(f3 * (1-f2));
   F_dummy9->matrix_multiply(minus_p_4, p_1);
   F_dummy9->multiply_by(f2 * (1-f3));
   F_dummy9->multiply_by(complex<double> (-1,0));
   F_dummy10->matrix_add(F_dummy8, F_dummy9);
   F_dummy10->multiply_by(pow(_sin_squared_theta_W_,2));
   
   F_dummy11->matrix_add(F_dummy7, F_dummy10);
   
   complex<double> comp_F0 = F_dummy11->get_A0();
   complex_three_vector* comp_F = F_dummy11->get_A();
   
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
   delete F_dummy10;
   delete F_dummy11;
   delete p_4;
   delete minus_p_4;
   delete p_1;
   delete minus_p_1;  
   delete A;
   delete G_R;
   delete G_L;
}

void nu_e_collision::F_LR_F_RL(double* F0, three_vector* F, density* dens, bool neutrino, int q2, double E2, int q3, double E3, int p4, double p4_energy, int count_min, int count_max){
   complex_three_vector* A = new complex_three_vector();
   A->set_value(2, complex<double> (0.5,0));
   matrix* G_L = new matrix(complex<double> (_sin_squared_theta_W_,0),A);
   
   matrix* G_R = new matrix(true);
   G_R->multiply_by(_sin_squared_theta_W_);
   
   matrix* p_1 = new matrix();
   matrix* minus_p_1 = new matrix();
   matrix* p_4 = new matrix();
   matrix* minus_p_4 = new matrix(true);
   
   p_1->convert_p_to_matrix(dens, neutrino, p1);
   minus_p_1->convert_p_to_identity_minus_matrix(dens, neutrino, p1);
   
   double f2 = 1 / (exp(E2/Tcm)+1);
   double f3 = 1 / (exp(E3/Tcm)+1);
   
    /*
   if(p4 == -1){
       p_4->convert_p4_to_interpolated_matrix(dens, neutrino, p4_energy);
       minus_p_4->convert_p4_to_identity_minus_interpolated_matrix(dens, neutrino, p4_energy);
   }
   else if(p4 == -2){
       p_4->convert_p4_to_interpolated_matrix(dens, neutrino, p4_energy);
       minus_p_4->convert_p4_to_identity_minus_interpolated_matrix(dens, neutrino, p4_energy);
   }
   else{
       p_4->convert_p_to_matrix(dens, neutrino, p4);
       minus_p_4->convert_p_to_identity_minus_matrix(dens, neutrino, p4);
   }*/
   
/*
   F_dummy1 = G_L * rho_4
   F_dummy2 = F_dummy1 * (1-rho_1)
   F_dummy3 = rho_4 * G_L
   F_dummy4 = F_dummy3 * (1-rho_1)
   
   F_dummy5 = G_L * (1-rho_4)
   F_dummy6 = F_dumy5 * rho_1
   F_dummy7 = (1-rho_4) * G_L
   F_dummy8 = F_dummy7 * rho_1
   
   F_dummy9 = F_dummy2 + F_dummy4 => F_dummy9=F_LR
   F_dummy10 = F_dummy6 + F_dummy8 => F_dummy10=F_RL
   F_dummy11 = F_dummy9 - F_dummy10 => F_dummy11 = F_LR + F_RL
   */

   matrix* F_dummy1 = new matrix();
   matrix* F_dummy2 = new matrix();
   matrix* F_dummy3 = new matrix();
   matrix* F_dummy4 = new matrix();
   matrix* F_dummy5 = new matrix();
   matrix* F_dummy6 = new matrix();
   matrix* F_dummy7 = new matrix();
   matrix* F_dummy8 = new matrix();
   matrix* F_dummy9 = new matrix();
   matrix* F_dummy10 = new matrix();
   matrix* F_dummy11 = new matrix();
   
   F_dummy1->matrix_multiply(G_L, p_4);
   F_dummy2->matrix_multiply(F_dummy1, minus_p_1);
   F_dummy3->matrix_multiply(p_4, G_L);
   F_dummy4->matrix_multiply(F_dummy3, minus_p_1);
   
   F_dummy5->matrix_multiply(G_L, minus_p_4);
   F_dummy6->matrix_multiply(F_dummy5, p_1);
   F_dummy7->matrix_multiply(minus_p_4, G_L);
   F_dummy8->matrix_multiply(F_dummy7, p_1);
   
   F_dummy9->matrix_add(F_dummy2, F_dummy4);
   F_dummy9->multiply_by(f3 * (1-f2));
   F_dummy10->matrix_add(F_dummy6, F_dummy8);
   F_dummy10->multiply_by(f2 * (1-f3));
   F_dummy10->multiply_by(complex<double> (-1,0));
   F_dummy11->matrix_add(F_dummy9, F_dummy10);
   F_dummy11->multiply_by(_sin_squared_theta_W_);
   
   complex<double> comp_F0 = F_dummy11->get_A0();
   complex_three_vector* comp_F = F_dummy11->get_A();
   
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
   delete F_dummy10;
   delete F_dummy11;
   delete p_4;
   delete minus_p_4;
   delete p_1;
   delete minus_p_1;  
   delete A;
   delete G_R;
   delete G_L;
}

void nu_e_collision::all_F_for_p1(density* dens, bool neutrino){
   double F0 = 0;
   three_vector* Fxyz = new three_vector();
   
   int p4_index = 0;
   double p4_energy = 0;
   double E3 = 0;
   double E2 = 0;
   for(int q3=0; q3<q3_vals_R2->get_len(); q3++){
       E3 = sqrt(pow(q3_vals_R2->get_value(q3),2) + me_squared);
       for(int q2=0; q2<q2_vals_R2[q3]->get_len(); q2++){
           E2 = sqrt(pow(q2_vals_R2[q3]->get_value(q2),2) + me_squared);
           p4_energy = p1_energy + E2 - E3;
           if(q2==0){
               //this q2 corresponds to q2_min so p4_energy won't be in eps so p4_index doesn't mean anything and instead is just an indicator of this case (p4_energy is minimized)
               p4_index = -1;
           }
           else if(q2==q2_vals_R2[q3]->get_len()-1){
               //this q2 corresponds to q2_max so p4_energy won't be in eps so p4_index doesn't mean anything and instead is just an indicator of this case (p4 energy is maximized)
               p4_index = -2;
           }
           else{
               //because count_min should give index of p4 energy that corresponds to q2_vals_R2[1]
               //note that by its construction count_min >= 1 so p4_index>=0
               p4_index = count_min_vals_R2[q3]+q2-1;
           }
           
           if(q2>=eps->get_len()+1){
               std::cout << "going to get a segmentation fault" << std::endl;
           }
           
           F_LR_F_RL(&F0, Fxyz, dens, neutrino, q2, E2, q3, E3, p4_index, p4_energy, count_min_vals_R2[q3], count_max_vals_R2[q3]);
               
           R2_F_LR_RL_values[0][q2][q3] = F0;
           R2_F_LR_RL_values[1][q2][q3] = Fxyz->get_value(0);
           R2_F_LR_RL_values[2][q2][q3] = Fxyz->get_value(1);
           R2_F_LR_RL_values[3][q2][q3] = Fxyz->get_value(2);
           
           F_LL_F_RR(&F0, Fxyz, dens, neutrino, q2, E2, q3, E3, p4_index, p4_energy, count_min_vals_R2[q3], count_max_vals_R2[q3]);
           
           R2_F_LL_RR_values[0][q2][q3] = F0;
           R2_F_LL_RR_values[1][q2][q3] = Fxyz->get_value(0);
           R2_F_LL_RR_values[2][q2][q3] = Fxyz->get_value(1);
           R2_F_LL_RR_values[3][q2][q3] = Fxyz->get_value(2);
       }
   }
   for(int q2=0; q2<q2_vals_R1->get_len(); q2++){
       E2 = sqrt(pow(q2_vals_R1->get_value(q2),2) + me_squared);
       for(int q3=0; q3<q3_vals_R1[q2]->get_len(); q3++){
           E3 = sqrt(pow(q3_vals_R1[q2]->get_value(q3),2) + me_squared);
           p4_energy = p1_energy + E2 - E3;
           if(q3==0){
               //this q3 corresponds to q3_min so p4_energy won't be in eps so p4_index doesn't mean anything and instead is just an indicator of this case (p4 energy maximized)
               p4_index = -2;
           }
           else if(q3==q3_vals_R1[q2]->get_len()-1){
               //this q3 corresponds to q3_max so p4_energy won't be in eps so p4_index doesn't mean anything and instead is just an indicator of this case (p4 energy minimized)
               p4_index = -1;
           }
           else{
               //because count_min should give index of p4 energy that corresponds to q3_vals_R1[1]
               //note that by its construction count_min >= 1 so p4_index>=0
               p4_index = count_min_vals_R1[q2]+q3-1;
           }
           
           if(q3>=eps->get_len()+1){
               std::cout << "going to get a segmentation fault, trying to plug q3 index of " << q3 << std::endl;
           }
           
           F_LR_F_RL(&F0, Fxyz, dens, neutrino, q2, E2, q3, E3, p4_index, p4_energy, count_min_vals_R1[q2], count_max_vals_R1[q2]);
               
           R1_F_LR_RL_values[0][q2][q3] = F0;
           R1_F_LR_RL_values[1][q2][q3] = Fxyz->get_value(0);
           R1_F_LR_RL_values[2][q2][q3] = Fxyz->get_value(1);
           R1_F_LR_RL_values[3][q2][q3] = Fxyz->get_value(2);
           
           F_LL_F_RR(&F0, Fxyz, dens, neutrino, q2, E2, q3, E3, p4_index, p4_energy, count_min_vals_R1[q2], count_max_vals_R1[q2]);
           
           R1_F_LL_RR_values[0][q2][q3] = F0;
           R1_F_LL_RR_values[1][q2][q3] = Fxyz->get_value(0);
           R1_F_LL_RR_values[2][q2][q3] = Fxyz->get_value(1);
           R1_F_LL_RR_values[3][q2][q3] = Fxyz->get_value(2);
       }
   }
   
   delete Fxyz;    
}

double nu_e_collision::M_11(int which, double q2_momentum, double q3_momentum, double E2, double E3){
   double b = 0;
   double a = 0;
   double C1 = pow(p1_energy + E2,2) - me_squared;
   
   if(which==1){
       a = p1_energy + E2 - E3 - q3_momentum;
       b = p1_energy + E2 - E3 + q3_momentum;
   }
   else if(which==2){
       a = p1_energy - q2_momentum;
       b = p1_energy + q2_momentum;
   }
   else if(which==3){
       a = E3 + q3_momentum - p1_energy - E2;
       b = p1_energy + q2_momentum;
   }
   else{
       a = q2_momentum - p1_energy;
       b = p1_energy + E2 - E3 + q3_momentum;
   }
   
   return 0.5 * (C1 * b - 1./3 * pow(b,3)) - 0.5 * (C1 * a - 1./3 * pow(a,3));
}

double nu_e_collision::M_12(int which, double q2_momentum, double q3_momentum, double E2, double E3){
   double b = 0;
   double a = 0;
   double C1 = pow(p1_energy + E2,2) - me_squared;
   
   if(which==1){
       a = p1_energy + E2 - E3 - q3_momentum;
       b = p1_energy + E2 - E3 + q3_momentum;
   }
   else if(which==2){
       a = p1_energy - q2_momentum;
       b = p1_energy + q2_momentum;
   }
   else if(which==3){
       a = E3 + q3_momentum - p1_energy - E2;
       b = p1_energy + q2_momentum;
   }
   else{
       a = q2_momentum - p1_energy;
       b = p1_energy + E2 - E3 + q3_momentum;
   }
   
   return 0.25 * (pow(C1,2) * b - 2./3 * C1 * pow(b,3) + 1./5 * pow(b,5)) - 0.25 * (pow(C1,2) * a - 2./3 * C1 * pow(a,3) + 1./5 * pow(a,5));
   
}

double nu_e_collision::M_21(int which, double q2_momentum, double q3_momentum, double E2, double E3){
   double b = 0;
   double a = 0;
   double C2 = pow(p1_energy - E3,2) - me_squared;
   
   if(which==1){
       a = p1_energy + E2 - E3 - q2_momentum;
       b = p1_energy + E2 - E3 + q2_momentum;
   }
   else if(which==2){
       a = p1_energy - q3_momentum;
       b = p1_energy + q3_momentum;
   }
   else if(which==3){
       a = E3 - p1_energy - E2 + q2_momentum;
       b = p1_energy + q3_momentum;
   }
   else{
       a = q3_momentum - p1_energy;
       b = p1_energy + E2 - E3 + q2_momentum;
   }
   
   return 0.5 * (1./3 * pow(b,3) - C2 * b) - 0.5 * (1./3 * pow(a,3) - C2 * a);
   
}

double nu_e_collision::M_22(int which, double q2_momentum, double q3_momentum, double E2, double E3){
   double b = 0;
   double a = 0;
   double C2 = pow(p1_energy - E3,2) - me_squared;
   
   if(which==1){
       a = p1_energy + E2 - E3 - q2_momentum;
       b = p1_energy + E2 - E3 + q2_momentum;
   }
   else if(which==2){
       a = p1_energy - q3_momentum;
       b = p1_energy + q3_momentum;
   }
   else if(which==3){
       a = E3 - p1_energy - E2 + q2_momentum;
       b = p1_energy + q3_momentum;
   }
   else{
       a = q3_momentum - p1_energy;
       b = p1_energy + E2 - E3 + q2_momentum;
   }
   
   return 0.25 * (1./5 * pow(b,5) - 2./3 * C2 * pow(b,3) + pow(C2,2) * b) - 0.25 * (1./5 * pow(a,5) - 2./3 * C2 * pow(a,3) + pow(C2,2) * a);
}

double nu_e_collision::R2_inner_integral(int which_term, int q3){
   double q3_momentum = q3_vals_R2->get_value(q3);
   double E3 = sqrt(pow(q3_momentum,2) + me_squared);
   
   int term = 0;
   //case 1
   if(p1_me < (sqrt(5)-1)/4.){
       //case 1a
       if(q3_momentum < q_cut_3){
           for(int q2=0; q2<q2_vals_R2[q3]->get_len(); q2++){
               double q2_momentum = q2_vals_R2[q3]->get_value(q2);
               double E2 = sqrt(pow(q2_momentum,2) + me_squared);
               
               //case 1ai
               if(q2_momentum < q3_momentum){
                   term = 1;
               }
               //case 1aii
               else if(q2_momentum < q_trans_2_R2->get_value(q3)){
                   term = 2;
               }
               //case 1aiii
               else{
                   term = 3;
               } 
               
               inner_vals_R2[q3]->set_value(q2, q2_momentum / E2 * (4 * R2_F_LL_RR_values[which_term][q2][q3] * M_22(term, q2_momentum, q3_momentum, E2, E3) + 4 * me_squared * R2_F_LR_RL_values[which_term][q2][q3] * M_21(term, q2_momentum, q3_momentum, E2, E3)));
               
           }
       }
       //case 1b
       else if(q3_momentum < q_cut_2_R2){
           for(int q2=0; q2<q2_vals_R2[q3]->get_len(); q2++){
               double q2_momentum = q2_vals_R2[q3]->get_value(q2);
               double E2 = sqrt(pow(q2_momentum,2) + me_squared);
               
               //case 1bi
               if(q2_momentum < q_trans_2_R2->get_value(q3)){
                   term = 1;
               }
               //case 1bii
               else if(q2_momentum < q3_momentum){
                   term = 4;
               }
               //case 1biii
               else{
                   term = 3;
               }
               
               inner_vals_R2[q3]->set_value(q2, q2_momentum / E2 * (4 * R2_F_LL_RR_values[which_term][q2][q3] * M_22(term, q2_momentum, q3_momentum, E2, E3) + 4 * me_squared * R2_F_LR_RL_values[which_term][q2][q3] * M_21(term, q2_momentum, q3_momentum, E2, E3)));
           }
       }
       //case 1c
       else if(q3_momentum < q_cut_1_R2){
           for(int q2=0; q2<q2_vals_R2[q3]->get_len(); q2++){
               double q2_momentum = q2_vals_R2[q3]->get_value(q2);
               double E2 = sqrt(pow(q2_momentum,2) + me_squared);
               
               //case 1ci
               if(q2_momentum < q3_momentum){
                   term = 4;
               }
               //case 1cii
               else{
                   term = 3;
               }
               
               inner_vals_R2[q3]->set_value(q2, q2_momentum / E2 * (4 * R2_F_LL_RR_values[which_term][q2][q3] * M_22(term, q2_momentum, q3_momentum, E2, E3) + 4 * me_squared * R2_F_LR_RL_values[which_term][q2][q3] * M_21(term, q2_momentum, q3_momentum, E2, E3)));
           }
       }
       //case 1d
       else{
           for(int q2=0; q2<q2_vals_R2[q3]->get_len(); q2++){
               double q2_momentum = q2_vals_R2[q3]->get_value(q2);
               double E2 = sqrt(pow(q2_momentum,2) + me_squared);
               
               //case 1di
               if(q2_momentum < q3_momentum){
                   term = 4;
               }
               //case 1dii
               else{
                   term = 3;
               }
               
               inner_vals_R2[q3]->set_value(q2, q2_momentum / E2 * (4 * R2_F_LL_RR_values[which_term][q2][q3] * M_22(term, q2_momentum, q3_momentum, E2, E3) + 4 * me_squared * R2_F_LR_RL_values[which_term][q2][q3] * M_21(term, q2_momentum, q3_momentum, E2, E3)));
           }
       }
   }
   
   //case 2
   else if(p1_me < 1./(2 * sqrt(2))){
       //case 2a
       if(q3_momentum<q_cut_3){
           for(int q2=0; q2<q2_vals_R2[q3]->get_len(); q2++){
               double q2_momentum = q2_vals_R2[q3]->get_value(q2);
               double E2 = sqrt(pow(q2_momentum,2) + me_squared);
               
               //case 2ai
               if(q2_momentum < q3_momentum){
                   term = 1;
               }
               //case 2aii
               else if(q2_momentum < q_trans_2_R2->get_value(q3)){
                   term = 2;
               }
               //case 2aiii
               else{
                   term = 3;
               }
               
               inner_vals_R2[q3]->set_value(q2, q2_momentum / E2 * (4 * R2_F_LL_RR_values[which_term][q2][q3] * M_22(term, q2_momentum, q3_momentum, E2, E3) + 4 * me_squared * R2_F_LR_RL_values[which_term][q2][q3] * M_21(term, q2_momentum, q3_momentum, E2, E3)));
           }
       }
       //case 2b
       else if(q3_momentum<q_cut_1_R2){
           for(int q2=0; q2<q2_vals_R2[q3]->get_len(); q2++){
               double q2_momentum = q2_vals_R2[q3]->get_value(q2);
               double E2 = sqrt(pow(q2_momentum,2) + me_squared);
               
               //case 2bi
               if(q2_momentum < q_trans_2_R2->get_value(q3)){
                   term = 1;
               }
               //case 2bii
               else if(q2_momentum < q3_momentum){
                   term = 4;
               }
               //case 2biii
               else{
                   term = 3;
               }
               
               inner_vals_R2[q3]->set_value(q2, q2_momentum / E2 * (4 * R2_F_LL_RR_values[which_term][q2][q3] * M_22(term, q2_momentum, q3_momentum, E2, E3) + 4 * me_squared * R2_F_LR_RL_values[which_term][q2][q3] * M_21(term, q2_momentum, q3_momentum, E2, E3)));
           }
       }
       //case 2c
       else if(q3_momentum<q_cut_2_R2){
           for(int q2=0; q2<q2_vals_R2[q3]->get_len(); q2++){
               double q2_momentum = q2_vals_R2[q3]->get_value(q2);
               double E2 = sqrt(pow(q2_momentum,2) + me_squared);
               
               //case 2ci
               if(q2_momentum < q_trans_2_R2->get_value(q3)){
                   term = 1;
               }
               //case 2cii
               else if(q2_momentum < q3_momentum){
                   term = 4;
               }
               //case 2ciii
               else{
                   term = 3;
               }
               
               inner_vals_R2[q3]->set_value(q2, q2_momentum / E2 * (4 * R2_F_LL_RR_values[which_term][q2][q3] * M_22(term, q2_momentum, q3_momentum, E2, E3) + 4 * me_squared * R2_F_LR_RL_values[which_term][q2][q3] * M_21(term, q2_momentum, q3_momentum, E2, E3)));
           }
       }
       //case 2d
       else{
           for(int q2=0; q2<q2_vals_R2[q3]->get_len(); q2++){
               double q2_momentum = q2_vals_R2[q3]->get_value(q2);
               double E2 = sqrt(pow(q2_momentum,2) + me_squared);
               
               //case 2di
               if(q2_momentum < q3_momentum){
                   term = 4;
               }
               //case 2dii
               else{
                   term = 3;
               }
               
               inner_vals_R2[q3]->set_value(q2, q2_momentum / E2 * (4 * R2_F_LL_RR_values[which_term][q2][q3] * M_22(term, q2_momentum, q3_momentum, E2, E3) + 4 * me_squared * R2_F_LR_RL_values[which_term][q2][q3] * M_21(term, q2_momentum, q3_momentum, E2, E3)));
           }
       }
       
   }
   
   //case 3
   else if(p1_me < 0.5){
       //case 3a
       if(q3_momentum<q_cut_1_R2){
           for(int q2=0; q2<q2_vals_R2[q3]->get_len(); q2++){
               double q2_momentum = q2_vals_R2[q3]->get_value(q2);
               double E2 = sqrt(pow(q2_momentum,2) + me_squared);
               
               //case 3ai
               if(q2_momentum < q3_momentum){
                   term = 1;
               }
               //case 3aii
               else if(q2_momentum < q_trans_2_R2->get_value(q3)){
                   term = 2;
               }
               //case 3aiii
               else{
                   term = 3;
               }
               
               inner_vals_R2[q3]->set_value(q2, q2_momentum / E2 * (4 * R2_F_LL_RR_values[which_term][q2][q3] * M_22(term, q2_momentum, q3_momentum, E2, E3) + 4 * me_squared * R2_F_LR_RL_values[which_term][q2][q3] * M_21(term, q2_momentum, q3_momentum, E2, E3)));
           }
       }
       //case 3b
       else if(q3_momentum<q_cut_3){
           for(int q2=0; q2<q2_vals_R2[q3]->get_len(); q2++){
               double q2_momentum = q2_vals_R2[q3]->get_value(q2);
               double E2 = sqrt(pow(q2_momentum,2) + me_squared);
               
               //case 3bi
               if(q2_momentum < q3_momentum){
                   term = 1;
               }
               //case 3bii
               else if(q2_momentum < q_trans_2_R2->get_value(q3)){
                   term = 2;
               }
               //case 3biii
               else{
                   term = 3;
               }
               
               inner_vals_R2[q3]->set_value(q2, q2_momentum / E2 * (4 * R2_F_LL_RR_values[which_term][q2][q3] * M_22(term, q2_momentum, q3_momentum, E2, E3) + 4 * me_squared * R2_F_LR_RL_values[which_term][q2][q3] * M_21(term, q2_momentum, q3_momentum, E2, E3)));
           }
       }
       //case 3c
       else if(q3_momentum<q_cut_2_R2){
           for(int q2=0; q2<q2_vals_R2[q3]->get_len(); q2++){
               double q2_momentum = q2_vals_R2[q3]->get_value(q2);
               double E2 = sqrt(pow(q2_momentum,2) + me_squared);
               
               //case 3ci
               if(q2_momentum < q_trans_2_R2->get_value(q3)){
                   term = 1;
               }
               //case 3cii
               else if(q2_momentum < q3_momentum){
                   term = 4;
               }
               //case 3ciii
               else{
                   term = 3;
               }
               
               inner_vals_R2[q3]->set_value(q2, q2_momentum / E2 * (4 * R2_F_LL_RR_values[which_term][q2][q3] * M_22(term, q2_momentum, q3_momentum, E2, E3) + 4 * me_squared * R2_F_LR_RL_values[which_term][q2][q3] * M_21(term, q2_momentum, q3_momentum, E2, E3)));
           }
       }
       //case 3d
       else{
           for(int q2=0; q2<q2_vals_R2[q3]->get_len(); q2++){
               double q2_momentum = q2_vals_R2[q3]->get_value(q2);
               double E2 = sqrt(pow(q2_momentum,2) + me_squared);
               
               //case 3di
               if(q2_momentum < q3_momentum){
                   term = 4;
               }
               //case 3dii
               else{
                   term = 3;
               }
               
               inner_vals_R2[q3]->set_value(q2, q2_momentum / E2 * (4 * R2_F_LL_RR_values[which_term][q2][q3] * M_22(term, q2_momentum, q3_momentum, E2, E3) + 4 * me_squared * R2_F_LR_RL_values[which_term][q2][q3] * M_21(term, q2_momentum, q3_momentum, E2, E3)));
           }
       }        
   }
   
   //case 4
   else{
        //case 4a
       if(q3_momentum<q_cut_1_R2){
           for(int q2=0; q2<q2_vals_R2[q3]->get_len(); q2++){
               double q2_momentum = q2_vals_R2[q3]->get_value(q2);
               double E2 = sqrt(pow(q2_momentum,2) + me_squared);
               
               //case 4ai
               if(q2_momentum < q3_momentum){
                   term = 1;
               }
               //case 4aii
               else{
                   term = 2;
               }
               
               inner_vals_R2[q3]->set_value(q2, q2_momentum / E2 * (4 * R2_F_LL_RR_values[which_term][q2][q3] * M_22(term, q2_momentum, q3_momentum, E2, E3) + 4 * me_squared * R2_F_LR_RL_values[which_term][q2][q3] * M_21(term, q2_momentum, q3_momentum, E2, E3)));
           } 
       }
       //case 4b
       else if(q3_momentum<q_cut_3){
           for(int q2=0; q2<q2_vals_R2[q3]->get_len(); q2++){
               double q2_momentum = q2_vals_R2[q3]->get_value(q2);
               double E2 = sqrt(pow(q2_momentum,2) + me_squared);
               
               //case 4bi
               if(q2_momentum < q3_momentum){
                   term = 1;
               }
               //case 4bii
               else if(q2_momentum < q_trans_2_R2->get_value(q3)){
                   term = 2;
               }
               //case 4biii
               else{
                   term = 3;
               }
               
               inner_vals_R2[q3]->set_value(q2, q2_momentum / E2 * (4 * R2_F_LL_RR_values[which_term][q2][q3] * M_22(term, q2_momentum, q3_momentum, E2, E3) + 4 * me_squared * R2_F_LR_RL_values[which_term][q2][q3] * M_21(term, q2_momentum, q3_momentum, E2, E3)));
           }
       }
       //case 4c
       else if(q3_momentum<q_cut_2_R2){
           for(int q2=0; q2<q2_vals_R2[q3]->get_len(); q2++){
               double q2_momentum = q2_vals_R2[q3]->get_value(q2);
               double E2 = sqrt(pow(q2_momentum,2) + me_squared);
               
               //case 4ci
               if(q2_momentum < q_trans_2_R2->get_value(q3)){
                   term = 1;
               }
               //case 4cii
               else if(q2_momentum < q3_momentum){
                   term = 4;
               }
               //case 4ciii
               else{
                   term = 3;
               }
               
               inner_vals_R2[q3]->set_value(q2, q2_momentum / E2 * (4 * R2_F_LL_RR_values[which_term][q2][q3] * M_22(term, q2_momentum, q3_momentum, E2, E3) + 4 * me_squared * R2_F_LR_RL_values[which_term][q2][q3] * M_21(term, q2_momentum, q3_momentum, E2, E3)));
           }
       }
       //case 4d
       else{
           for(int q2=0; q2<q2_vals_R2[q3]->get_len(); q2++){
               double q2_momentum = q2_vals_R2[q3]->get_value(q2);
               double E2 = sqrt(pow(q2_momentum,2) + me_squared);
               
               //case 4di
               if(q2_momentum < q3_momentum){
                   term = 4;
               }
               //case 4dii
               else{
                   term = 3;
               }
               
               inner_vals_R2[q3]->set_value(q2, q2_momentum / E2 * (4 * R2_F_LL_RR_values[which_term][q2][q3] * M_22(term, q2_momentum, q3_momentum, E2, E3) + 4 * me_squared * R2_F_LR_RL_values[which_term][q2][q3] * M_21(term, q2_momentum, q3_momentum, E2, E3)));
           }
           
       }  
       
                   
   }
   double result = q2_vals_R2[q3]->integrate(inner_vals_R2[q3]);
   return result;
}

void nu_e_collision::R2_whole_integral(double* results){
   double q3_momentum = 0;
   double E3 = 0;
   for(int i=0; i<4; i++){
       for(int q3=0; q3<q3_vals_R2->get_len(); q3++){
           q3_momentum = q3_vals_R2->get_value(q3);
           E3 = sqrt(pow(q3_momentum,2) + me_squared);
           outer_vals_R2->set_value(q3, q3_momentum / E3 * R2_inner_integral(i, q3));     
       }
       results[i] = q3_vals_R2->integrate(outer_vals_R2);
       results[i] *= pow(Tcm,5) / (pow(2,4) * pow(2*_PI_,3) * pow(p1_energy,2));
   }                 
}


double nu_e_collision::R1_inner_integral(int which_term, int q2){
   double q2_momentum = q2_vals_R1->get_value(q2);
   double E2 = sqrt(pow(q2_momentum,2) + me_squared);
   
   int term = 0;
   //case 1
   if(p1_me < 0.5){
       //case 1a
       if(q2_momentum < q_cut_3){
           for(int q3=0; q3<q3_vals_R1[q2]->get_len(); q3++){
               double q3_momentum = q3_vals_R1[q2]->get_value(q3);
               double E3 = sqrt(pow(q3_momentum,2) + me_squared);
               
               //case 1ai
               if(q3_momentum < q2_momentum){
                   term = 1;
               }
               //case 1aii
               else if(q3_momentum < q_trans_2_R1->get_value(q2)){
                   term = 2;
               }
               //case 1aiii
               else{
                   term = 3;
               } 
               
               inner_vals_R1[q2]->set_value(q3, q3_momentum / E3 * (4 * R1_F_LL_RR_values[which_term][q2][q3] * M_12(term, q2_momentum, q3_momentum, E2, E3) - 4 * me_squared * R1_F_LR_RL_values[which_term][q2][q3] * M_11(term, q2_momentum, q3_momentum, E2, E3)));
               
           }
       }
       //case 1b
       else if(q2_momentum < q_cut_1_R1){
           for(int q3=0; q3<q3_vals_R1[q2]->get_len(); q3++){
               double q3_momentum = q3_vals_R1[q2]->get_value(q3);
               double E3 = sqrt(pow(q3_momentum,2) + me_squared);
               
               //case 1bi
               if(q3_momentum < q_trans_2_R1->get_value(q2)){
                   term = 1;
               }
               //case 1bii
               else if(q3_momentum < q2_momentum){
                   term = 4;
               }
               //case 1biii
               else{
                   term = 3;
               } 
               
               inner_vals_R1[q2]->set_value(q3, q3_momentum / E3 * (4 * R1_F_LL_RR_values[which_term][q2][q3] * M_12(term, q2_momentum, q3_momentum, E2, E3) - 4 * me_squared * R1_F_LR_RL_values[which_term][q2][q3] * M_11(term, q2_momentum, q3_momentum, E2, E3)));
               
           }
           
       }
       else{
           for(int q3=0; q3<q3_vals_R1[q2]->get_len(); q3++){
               double q3_momentum = q3_vals_R1[q2]->get_value(q3);
               double E3 = sqrt(pow(q3_momentum,2) + me_squared);
               
               //case 1ci
               if(q3_momentum < q2_momentum){
                   term = 4;
               }
               //case 1cii
               else{
                   term = 3;
               } 
               
               inner_vals_R1[q2]->set_value(q3, q3_momentum / E3 * (4 * R1_F_LL_RR_values[which_term][q2][q3] * M_12(term, q2_momentum, q3_momentum, E2, E3) - 4 * me_squared * R1_F_LR_RL_values[which_term][q2][q3] * M_11(term, q2_momentum, q3_momentum, E2, E3)));
               
           }
       }
   }
   else{
       //case 2a
       if(q2_momentum < q_cut_3){
           for(int q3=0; q3<q3_vals_R1[q2]->get_len(); q3++){
               double q3_momentum = q3_vals_R1[q2]->get_value(q3);
               double E3 = sqrt(pow(q3_momentum,2) + me_squared);
               
               //case 2ai
               if(q3_momentum < q2_momentum){
                   term = 1;
               }
               //case 2aii
               else if(q3_momentum < q_trans_2_R1->get_value(q2)){
                   term = 2;
               }
               //case 2aiii
               else{
                   term = 3;
               } 
               
               inner_vals_R1[q2]->set_value(q3, q3_momentum / E3 * (4 * R1_F_LL_RR_values[which_term][q2][q3] * M_12(term, q2_momentum, q3_momentum, E2, E3) - 4 * me_squared * R1_F_LR_RL_values[which_term][q2][q3] * M_11(term, q2_momentum, q3_momentum, E2, E3)));
               
           }
           
       }
       //case 2b
       else{
           for(int q3=0; q3<q3_vals_R1[q2]->get_len(); q3++){
               double q3_momentum = q3_vals_R1[q2]->get_value(q3);
               double E3 = sqrt(pow(q3_momentum,2) + me_squared);
               
               //case 2bi
               if(q3_momentum < q_trans_2_R1->get_value(q2)){
                   term = 1;
               }
               //case 2bii
               else if(q3_momentum < q2_momentum){
                   term = 4;
               }
               //case 2biii
               else{
                   term = 3;
               } 
               
               inner_vals_R1[q2]->set_value(q3, q3_momentum / E3 * (4 * R1_F_LL_RR_values[which_term][q2][q3] * M_12(term, q2_momentum, q3_momentum, E2, E3) - 4 * me_squared * R1_F_LR_RL_values[which_term][q2][q3] * M_11(term, q2_momentum, q3_momentum, E2, E3)));
           }
       }
   }
   double result = q3_vals_R1[q2]->integrate(inner_vals_R1[q2]);
   
   return result; 
   
}

void nu_e_collision::R1_whole_integral(double* results){
   double q2_momentum = 0;
   double E2 = 0;
   for(int i=0; i<4; i++){
       for(int q2=0; q2<q2_vals_R1->get_len(); q2++){
           q2_momentum = q2_vals_R1->get_value(q2);
           E2 = sqrt(pow(q2_momentum,2) + me_squared);
           outer_vals_R1->set_value(q2, q2_momentum / E2 * R1_inner_integral(i, q2));     
       }
       results[i] = q2_vals_R1->integrate(outer_vals_R1);
       results[i] *= pow(Tcm,5) / (pow(2,4) * pow(2*_PI_,3) * pow(p1_energy,2));
   }
}

void nu_e_collision::whole_integral(density* dens, bool neutrino, double* results){
   if(p1==0){
       for(int i=0; i<4; i++){
           results[i]=0;
       }
   }
   
   else{
       all_F_for_p1(dens, neutrino);
       double* results1 = new double[4]();
       double* results2 = new double[4]();
       R1_whole_integral(results1);
       R2_whole_integral(results2);

       for(int i=0; i<4; i++){
           results[i] = results1[i] + results2[i];
       }

       delete[] results1;
       delete[] results2;
   }
}


nu_e_collision::~nu_e_collision(){
   for(int i=0; i<4; i++){
       for(int j=0; j<eps->get_len()+1; j++){
           delete[] R2_F_LR_RL_values[i][j];
           delete[] R2_F_LL_RR_values[i][j];
       }
       delete[] R2_F_LR_RL_values[i];
       delete[] R2_F_LL_RR_values[i];
   }
   delete[] R2_F_LR_RL_values;
   delete[] R2_F_LL_RR_values;
   
   for(int i=0; i<q3_vals_R2->get_len(); i++){
       delete q2_vals_R2[i];
       delete inner_vals_R2[i];
   }
   delete q2_vals_R2;
   delete inner_vals_R2;
   
   delete outer_vals_R2;
   delete q3_vals_R2;
   
   
   delete eps;
   delete q_trans_2_R2;
   delete q_lim_1_R2;   
}
