#include "arrays.hh"
#include "matrices.hh"


class nu_nu_collision_one
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
    /*
    std::ofstream& p4s1;
    std::ofstream& p4vals1;
    std::ofstream& p4s2;
    std::ofstream& p4vals2;*/
    
    nu_nu_collision_one(linspace_and_gl*, int);
    
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
    
    ~nu_nu_collision_one();
    
};

class nu_nu_collision_two
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
    
    nu_nu_collision_two(linspace_and_gl*, int);
    
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
    
    ~nu_nu_collision_two();
    
};

class nu_nu_collision_one_1
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
    
    nu_nu_collision_one_1(linspace_and_gl*, int);
    
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
    
    ~nu_nu_collision_one_1();
    
};