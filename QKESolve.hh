#include "ODESolve.hh"
#include "arrays.hh"
#include "QKE_methods.hh"

class QKE : public ODESolve<density>
{

    protected:
        double delta_m_squared;
        double cos_2theta;
        double sin_2theta;

        double eta_e;
        double eta_mu;

        three_vector_for_QKE* dummy_v_vac;

    public:
        QKE(dummy_vars* epsilon, double cos_2theta, double delta_m_squared, double eta_e=0., double eta_mu=0.);
        ~QKE();
        void f(double, density*, density*);
};