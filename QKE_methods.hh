#ifndef _QKE_METHODS_HH
#define _QKE_METHODS_HH
#include "arrays.hh"
#include "matrices.hh"

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
    density(int, dummy_vars*, double*);
    density(dummy_vars*, double, double);
    density(dummy_vars*, int, int);
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
    
    double interpolate_p0(bool, double);
    void interpolate_p0p(bool, double, three_vector*);
};

class nu_nu_collision
{
    protected:
    linspace_and_gl* eps;
    int p1;
    dep_vars* outer_vals;
    dep_vars** inner_vals;
    dummy_vars** p3_vals;
    double*** Fvv_values;
    double*** Fvvbar_values;
    int count;
    
    public:
    
    nu_nu_collision(linspace_and_gl*, int);
    
    void Fvvsc_components_term_1(density*, bool, int, int, double*, three_vector*);
    void Fvvsc_components_term_2(density*, bool, int,int, double*, three_vector*);
    void Fvvsc_components(density*, bool, int, int, double*, three_vector*);
    void Fvvsc_for_p1(density*, bool);
    void Fvvbarsc_components_term_1(density*, bool, int, int, double*, three_vector*);
    void Fvvbarsc_components_term_2(density*, bool, int,int, double*, three_vector*);
    void Fvvbarsc_components(density*, bool, int, int, double*, three_vector*);
    void Fvvbarsc_for_p1(density*, bool);    
    double J1(double, double, double);
    double J2(double, double);
    double J3(double, double, double);
    double K1(double, double);
    double K2(double, double, double);
    double K3(double, double, double);
    double interior_integral(int, int);
    void whole_integral(density*, bool, double*);
    
    ~nu_nu_collision();
    
};

class nu_e_collision
{
    protected:
    linspace_and_gl* eps;
    int p1;
    double Tcm;
    double p1_energy;
    double scaled_me;
    double me_squared;
    double p1_me;
    double q_cut_3;
    
    double q_cut_1_R2;
    double q_cut_2_R2;
    dep_vars* q_trans_2_R2;
    dep_vars* q_lim_1_R2;
    
    double q_cut_1_R1;
    dep_vars* q_trans_2_R1;
    dep_vars* q_lim_1_R1;
    
    dummy_vars* q3_vals_R2;
    dep_vars* outer_vals_R2;
    dummy_vars** q2_vals_R2;
    dep_vars** inner_vals_R2;
    
    dummy_vars* q2_vals_R1;
    dep_vars* outer_vals_R1;
    dummy_vars** q3_vals_R1;
    dep_vars** inner_vals_R1;
    
    double* p4_min_vals_R2;
    double* p4_max_vals_R2;
    int* count_min_vals_R2;
    int* count_max_vals_R2;
    
    double* p4_min_vals_R1;
    double* p4_max_vals_R1;
    int* count_min_vals_R1;
    int* count_max_vals_R1;
    
    double*** R2_F_LL_RR_values;
    double*** R2_F_LR_RL_values;
    
    double*** R1_F_LL_RR_values;
    double*** R1_F_LR_RL_values;
    
    
    public:
    nu_e_collision(linspace_and_gl*, int, double);
    
    void all_F_for_p1(density*, bool);
    void F_LL_F_RR(double*, three_vector*, density*, bool, int, double, int, double, int, double, int, int);
    void F_LR_F_RL(double*, three_vector*, density*, bool, int, double, int, double, int, double, int, int);
    
    double R2_inner_integral(int, int);
    void R2_whole_integral(double*);
    
    double R1_inner_integral(int, int);
    void R1_whole_integral(double*);
    
    void whole_integral(density*, bool, double*);
    
    double M_11(int, double, double, double, double);
    double M_12(int, double, double, double, double);
    double M_21(int, double, double, double, double);
    double M_22(int, double, double, double, double);
    
    ~nu_e_collision();
    
    
};

#endif