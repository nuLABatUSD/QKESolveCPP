#include "QKESolve.hh"

class QKESolveMPI : public QKE{
    private:
    int myid;
    int numprocs;
    
    public:
    QKESolveMPI(int, int, linspace_and_gl*, double, double, double, double);

    void f(double, density*, density*);
    void RKCash_Karp(double, density*, double, double*, density*, density*);
    bool step_accept(density*, density*, density*, double, double*);
    bool RKCK_step(double, density*, double, double*, density*, double*);
    bool ODEOneRun(double x0, density* y0, double dx0, int N_step, int dN, double x_final, double* x, density* y, double* dx, const std::string& file_name, bool verbose = false);

    bool run(int N_step, int dN, double x_final, const std::string& file_name, bool verbose = false);

};