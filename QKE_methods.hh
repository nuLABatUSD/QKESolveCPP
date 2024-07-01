#ifndef _QKE_METHODS_HH
#define _QKE_METHODS_HH
#include "arrays.hh"
#include <iostream>

using std::ostream;

class three_vector_for_QKE : public three_vector
{
    public:
    void v_vacuum(double, double, double);
    void v_density(dummy_vars*, density*);
    void v_thermal(dummy_vars*, density*);
};

class density : public dep_vars
{
    protected:
    int N_bins;
    dummy_vars* E;
    
    public:
    
    density(int, dummy_vars*);
    density(dummy_vars*, double, double);
    density(density*);
    ~density();
    
    double get_E_value(int);
    dummy_vars* get_E();
    double get_T();
    double get_Tcm();
    int num_bins();

    void set_T(double);
    double p0(int, bool);
    void p_vector(int, bool, three_vector*);
    void p0_p(int, bool, three_vector*);

    void number_density(double*);
    void print_csv(ostream&);
};

#endif