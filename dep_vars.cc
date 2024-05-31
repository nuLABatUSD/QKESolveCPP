#include <iostream>
#include "arrays.hh"
#include "constants.hh"
#include <cmath>
#include "energy_density_and_pressure.hh"

using std::cout;
using std::endl;
using std::ostream;

//dep_vars
dep_vars::dep_vars(int size)
{
    N = size;
    values = new double[N]();
}
    
dep_vars::dep_vars(double* copy_me, int size)
{
    N = size;
    values = new double[N]();
    for (int i = 0; i < N; i++)
        values[i] = copy_me[i];
}
    
dep_vars::dep_vars(dep_vars* copy_me)
{
    N = copy_me->length();
    values = new double[N]();
    for (int i = 0; i < N; i++)
        values[i] = copy_me->values[i];
}
    
dep_vars::~dep_vars()
{delete[] values;}

int dep_vars::length()
{return N;}

double dep_vars::get_value(int i)
{return values[i];}

void dep_vars::set_value(int i, double v)
{values[i] = v;}
    
void dep_vars::print_all()
{
    for (int i = 0; i < N; i++)
        cout << values[i] << endl;
}
    
void dep_vars::print(int N_top = 3, int N_bot = 1)
{
    if (N <= N_top + N_bot)
        print_all();
    else if (N_top < 0 || N_bot < 0)
        print_all();
    else
    {
        for (int i = 0; i < N_top; i++)
            cout << values[i] << endl;
        cout << "..." << endl;
        for (int i = 0; i < N_bot; i++)
            cout << values[N - N_bot + i] << endl;
    }
}

void dep_vars::multiply_by(double scalar)
{
    for (int i = 0; i < N; ++i) 
    {
        values[i] *= scalar;
    }
}

void dep_vars::copy(dep_vars* z)
{
    
    for (int i = 0; i < N; ++i) 
    {
        values[i] = z -> get_value(i);
    }
}

void dep_vars::add_to(double c, dep_vars* z)
{

    for (int i = 0; i < N; ++i) {
        //values[i] += c*z[i];
        values[i] += c * z -> get_value(i);
    }

}

// three_vector
three_vector::three_vector(int Nv):dep_vars(3)
{;}

three_vector::three_vector(double x, double y, double z):dep_vars(3)
{
    values[0] = x;
    values[1] = y;
    values[2] = z;
}

three_vector::three_vector(double* copy_me):dep_vars(copy_me, 3)
{;}

three_vector::three_vector(three_vector* copy_me):dep_vars(copy_me)
{;}

void three_vector::add(three_vector* A, three_vector*B){
    values[0] = A->get_value(0) + B->get_value(0);
    values[1] = A->get_value(1) + B->get_value(1);
    values[2] = A->get_value(2) + B->get_value(2);
}

double three_vector::dot_with(three_vector* B)
{
    double dot = 0;
    for(int i = 0; i < 3; i++)
        dot += values[i] * B->get_value(i);
    return dot;
}

double three_vector::magnitude_squared()
{
    return dot_with(this);
}

double three_vector::magnitude()
{
    double sum = 0;
    for(int i =0; i < 3; i++)
        sum += pow(this->get_value(i),2);
    return sqrt(sum);
}

void three_vector::set_cross_product(three_vector* A, three_vector* B)
{
    values[0] = A->get_value(1) * B->get_value(2) - A->get_value(2) * B->get_value(1);
    values[1] = A->get_value(2) * B->get_value(0) - A->get_value(0) * B->get_value(2);
    values[2] = A->get_value(0) * B->get_value(1) - A->get_value(1) * B->get_value(0);
}

void three_vector::v_vacuum(){
    values[0] = _delta_m_squared_ / 2 * _sin_2theta_;
    values[2] = _delta_m_squared_ / 2 * _cos_2theta_;
}

void three_vector::v_thermal(dummy_vars* q, density* d){
    
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
        values[i] *= -8*sqrt(2)*_GF_/(3*pow(_Z_boson_,2));
    }
    
    double energy_dens = 0;
    double pressure = 0;
    energy_and_pressure(_electron_mass_, d->get_T(), &energy_dens, &pressure);
    values[2] += -2 * sqrt(2) * pow(_W_boson_,-2) * _GF_ * (energy_dens+pressure);
    
    
}



void three_vector::v_density(dummy_vars* q, density* d){
    
    
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
        values[i] *= sqrt(2)*_GF_ / (2 * pow(_PI_,2));
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

dummy_vars* density::get_E(){
    return E;
}

double density::get_E_value(int i){
    return E->get_value(i);
}

double density::get_T(){
    return values[N-2];
}

int density::num_bins()
{return N_bins;}

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





