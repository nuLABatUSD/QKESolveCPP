#include "QKESolve.hh"
#include "QKE_methods.hh"
#include "arrays.hh"

class QKESolveMPI : public ODESolve<density>
{
    private:
    int myid;
    int numprocs;
    double delta_m_squared;
    double cos_2theta;
    double sin_2theta;

    linspace_and_gl* epsilon;
    three_vector_for_QKE* dummy_v_vac;
    nu_nu_collision** int_objects;
    QKE* just_h;
    
    bool find_nu_e;
    
    public:
    QKESolveMPI(int, int, linspace_and_gl*, double, double, double, double, const std::string&);
    void just_neutrino_collision();
    ~QKESolveMPI();

    void f(double, density*, density*);
    double first_derivative(double, density*, density*, double, double*);
    void RKCash_Karp(double, density*, double, double*, density*, density*);
    bool step_accept(density*, density*, density*, double, double*);
    bool RKCK_step(double, density*, double, double*, density*, double*);
    bool ODEOneRun(double x0, density* y0, double dx0, int N_step, int dN, double x_final, double* x, density* y, double* dx, const std::string& file_name, bool verbose = false);

    bool run(int N_step, int dN, double x_final, const std::string& file_name, bool verbose = false);

};