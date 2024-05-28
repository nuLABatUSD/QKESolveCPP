#include <iostream>
#include <cmath>
#include "arrays.hh"
#include "constants.hh"

using std::cout;
using std::endl;

three_vector* V_density();
//double V_vvbar();
//double VT_vvbar();


int main(){
    return 0;
}
/*
three_vector* V_vac(){
    new three_vector* vac;
    vac->set_value(0,_sin_2theta_);
    vac->set_value(1,0.);
    vac->set_value(2,_cos_2theta_);
    
    vac->multiply_by(_delta_m_squared_ / 2);
    
    return vac;
}
*/


three_vector V_density(dummy_vars* q, density* d){
    three_vector* v_dens = new three_vector();
    int N = q->get_len();
    double* dx = q->get_dx();
    
    three_vector* dummy1 = new three_vector();
    three_vector* dummy2 = new three_vector();
    for(int i=1; i<N; i++){
        d->p0_p(i, true, dummy1);
        d->p0_p(i, false, dummy2);

        v_dens->set_value(pow(q->get_val(i),2) * (dummy1->get_value(0) - dummy2->get_value(0)),0);
        v_dens->set_value(pow(q->get_val(i),2) * (dummy1->get_value(1) - dummy2->get_value(1)),1);
        v_dens->set_value(pow(q->get_val(i),2) * (dummy1->get_value(2) - dummy2->get_value(2)),2);
        
        v_dens->multiply_by(0.5 * dx[i]);
        
        
    }
    delete dummy1;
    delete dummy2;
    v_dens->multiply_by(sqrt(2)*_GF_ / (2* pow(_PI_,2)));
    return v_dens;
}
