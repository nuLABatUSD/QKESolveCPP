#include "QKE_methods.hh"
#include <iostream>
#include "constants.hh"
#include <cmath>
#include "energy_density_and_pressure.hh"
#include "gl_vals.hh"

three_vector_for_QKE::three_vector_for_QKE(double cos, double mass):three_vector(3){
    _delta_m_squared_ = mass;
    _cos_2theta_ = cos;
    _sin_2theta_ = sqrt(1-pow(_cos_2theta_,2));
    
}

void three_vector_for_QKE::v_vacuum(){
    values[0] = _delta_m_squared_ / 2 * _sin_2theta_;
    values[2] = _delta_m_squared_ / 2 * _cos_2theta_;
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
    E = eps;
    
}

density::density(dummy_vars* eps, double eta_nu, double eta_mu):dep_vars(8*eps->N+2)
{
    N_bins = eps->N;
    E = eps;

    double fnu = 0;
    double fmu = 0;
    double fnubar = 0;
    double fmubar = 0;
    
    for (int i=0; i<N_bins; i++){
        fnu = 1 / (exp(eps->values[i] - eta_nu)+1);
        fmu = 1 / (exp(eps->values[i] - eta_mu)+1);
        values[4*i] = fnu + fmu;
        values[4*i+3] =  (fnu - fmu)/(fnu+fmu);
       
       fnubar = 1 / (exp(eps->values[i] + eta_nu)+1);
       fmubar = 1 / (exp(eps->values[i] + eta_mu)+1);
       values[4*N_bins + 4*i] = fnubar + fmubar;
       values[4*N_bins + 4*i+3] = (fnubar - fmubar)/(fnubar+fmubar);
      
    }
    
}

density::density(density* copy_me):dep_vars(copy_me)
{
    
    N_bins = copy_me->num_bins();
    E = copy_me->get_E();

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
}

void density::print_csv(ostream& os)
{
    double nd[4];
    number_density(nd);

    for(int i = 0; i < 3; i++)
        os << nd[i] << ", ";
    os << nd[3];
}
    
    