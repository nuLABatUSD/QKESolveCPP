#include <iostream>
#include "constants.hh"
#include "ODEexample.hh"
#include "arrays.hh"

using std::cout;
using std::endl;

expo::expo() : ODESolve()
{
    y_values = new dep_vars(4);
}

void expo::f(double x, dep_vars* y, dep_vars* z)
{
    z->set_value(0, y->get_value(1));
    z->set_value(1, y->get_value(2));
    z->set_value(2, y->get_value(3));
    z->set_value(3, pow(_PI_, 4) * y->get_value(0));
    return;
}

spin::spin() : ODESolve()
{
    B = new three_vector(0.8, 0.0, 0.6);
    omega = _PI_;
    y_values = new three_vector(3);
}

void spin::f(double x, three_vector* y, three_vector* z)
{
    z->set_cross_product(B, y);
    z->multiply_by(omega);
}

void spin::print_state()
{
    ODESolve::print_state();
    cout << "| L | = " << y_values->magnitude() << endl;
}

QKE::QKE(int x, dummy_vars* E) : ODESolve()
{
    y_values = new density(x, E);
}
QKE::QKE(dummy_vars* E, double x, double y) : ODESolve()
{
    y_values = new density(E, x, y);
}

void QKE::f(double t, density* d1, density* d2)
{
    dummy_vars* E = d1->get_E();
    three_vector* dummy_v_vac = new three_vector();
    three_vector* dummy_v_dens = new three_vector();
    three_vector* dummy_v_therm = new three_vector();
    three_vector* dummy_v_1 = new three_vector();
    three_vector* dummy_v_2 = new three_vector();
    three_vector* dummy_v_3 = new three_vector();
    three_vector* vcrossp = new three_vector();
    
    three_vector* p = new three_vector();
    
    dummy_v_vac->v_vacuum();
    dummy_v_dens->v_density(E, d1);
    dummy_v_therm->v_thermal(E, d1);
    
    
    for (int i=1; i< E->get_len(); i++){
        dummy_v_1->copy(dummy_v_vac);
        dummy_v_1->multiply_by(1./E->get_value(i));

        dummy_v_2->copy(dummy_v_therm);
        dummy_v_2->multiply_by(E->get_value(i));
        
        dummy_v_3->add(dummy_v_vac, dummy_v_dens);
        dummy_v_3->add(dummy_v_1, dummy_v_therm);       

        
        d1->p_vector(i,true,p);
        
        vcrossp->set_cross_product(dummy_v_3,p);
        
        
        d2->set_value(4*i+1, vcrossp->get_value(0));
        d2->set_value(4*i+2, vcrossp->get_value(1));
        d2->set_value(4*i+3, vcrossp->get_value(2));
        
        
        d1->p_vector(i,false,p);
        vcrossp->set_cross_product(dummy_v_3,p);
        
        
        d2->set_value(4*(E->get_len())+4*i+1, vcrossp->get_value(0));
        d2->set_value(4*(E->get_len())+4*i+2, vcrossp->get_value(1));
        d2->set_value(4*(E->get_len())+4*i+3, vcrossp->get_value(2));
        
    }
        
    delete dummy_v_vac;
    delete dummy_v_dens;
    delete dummy_v_therm;
    delete dummy_v_1;
    delete dummy_v_2;
    delete dummy_v_3;
    delete vcrossp;
    delete p;
    
}

int main()
{    
    
    linspace_for_trap* et = new linspace_for_trap(0.,20, 201);
    double eta_e = 0.2;
    double eta_mu = -0.02;
    QKE* sim = new QKE(et, eta_e, eta_mu);
    density* den = new density(et, eta_e, eta_mu);
    sim->set_ics(0, den, 0.1);

    sim->print_state();

    sim->run(100,10,5.0,"output2.csv");

    sim->print_state();

    cout << "======================" << endl;
    return 1;
    
    /*
    expo* sim = new expo;
    double ics[] = {1.0, 0.0, - _PI_ * _PI_, 0.0};

    dep_vars* y0 = new dep_vars(ics, 4);
    sim->set_ics(0, y0, 0.1);

    sim->print_state();

    sim->run(100,10,5.0,"output2.csv");

    sim->print_state();

    cout << "======================" << endl;

    spin* sim2 = new spin;
    three_vector* L0 = new three_vector(0., 0., 1.);

    sim2->set_ics(0, L0, 0.1);
    sim2->print_state();

    sim2->run(100, 10, 5.0, "output3.csv");
    sim2->print_state();
    return 1;*/
}