#include <iostream>
#include "arrays.hh"
#include "constants.hh"
#include <cmath>

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
three_vector::three_vector():dep_vars(3)
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




//density

density::density(int num):dep_vars(8*num)
{
    N_bins = num;
}

density::density(linspace* eps, double eta_nu, double eta_mu):dep_vars(8*eps->N)
{
    N_bins = eps->N;

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
       values[4*N_bins + 4*i] = fnubar - fmubar;
       values[4*N_bins + 4*i+3] = (fnubar - fmubar)/(fnubar+fmubar);
      
    }
    
}

int density::num_bins()
{return N_bins;}

void density::p_vector(int t, bool neutrino, three_vector* p)
{
    if(neutrino==true){
        for(int i=0; i<3; i++){
            p->set_value(i, values[4*t+i+1]);
        }}
    else{
        for(int i=0; i<3; i++){
            p->set_value(i, values[N/2+4*t+1+i]);
        }}
}

void density::p0_p(int t, bool neutrino, three_vector* p)
{
    if(neutrino==true){
        for(int i=0; i<3; i++){
            p->set_value(i,values[4*t+i+1]);
        }
        p->multiply_by(values[t]);
        
    }
    
    else{
        for(int i=0; i<3; i++){
            p->set_value(i, values[N/2+4*t+i+1]);
        }
        p->multiply_by(values[t]);
    }
}

three_vector density::v_density_integral(dummy_vars* q, density* d){
    //NOTE CONSTANT TERM NOT INCLUDED!
    values[0]=0;
    values[1]=0;
    values[2]=0;
    
    three_vector* dummy1 = new three_vector();
    three_vector* dummy2 = new three_vector();
    for(int i=0; i<N; i++){
        d->p0_p(i, true, dummy1);
        d->p0_p(i, false, dummy2);

        values[0] += pow(q->get_val(i),2) * (dummy1->get_value(0) - dummy2->get_value(0));
        values[1] += pow(q->get_val(i),2) * (dummy1->get_value(1) - dummy2->get_value(1));
        values[2] += pow(q->get_val(i),2) * (dummy1->get_value(2) - dummy2->get_value(2));
        
        
    }
    delete dummy1;
    delete dummy2;
    return values;
}



