#include "QKESolveMPI.hh"
#include "mpi.h"

using std::abs;

/*
TO RUN:

mpic++ test2.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc -std=c++11 -o wed
mpiexec -n 4 wed
*/

//QKESolveMPI::QKESolveMPI(int rank, int numranks, linspace_and_gl* epss, double cos_2theta, double delta_m_squared, double eta_e=0., double eta_mu=0.) : QKE(epss, cos_2theta, delta_m_squared, eta_e, eta_mu){
QKESolveMPI::QKESolveMPI(int rank, int numranks, linspace_and_gl* e, double sin2theta, double deltamsquared, double eta_e, double eta_mu, double x0, double dx0, const std::string& dens_input) : ODESolve()
{
    myid = rank;
    numprocs = numranks;
    
    epsilon = new linspace_and_gl(e);
    sin_2theta = sin2theta;
    cos_2theta = sqrt(1 - pow(sin2theta, 2));
    delta_m_squared = deltamsquared;

    dummy_v_vac = new three_vector_for_QKE;
    dummy_v_vac->v_vacuum(delta_m_squared, cos_2theta, sin_2theta);
    
    int_objects = new nu_nu_collision*[epsilon->get_len()];
    for(int i=0; i<epsilon->get_len(); i++){
        int_objects[i] = new nu_nu_collision(epsilon, i);
    }
    
    just_h = new QKE(e, sin2theta, deltamsquared, eta_e, eta_mu);
    
    
    double* dens_vals = new double[8*epsilon->get_len()+2]();
    
    std::ifstream densfile;
    densfile.open(dens_input);
    if (!densfile.is_open()) {
        std::cout << "Error opening density input file" << std::endl;
    }
    
    std::string mystring;
    int i=0;
    while(densfile){
        std::getline(densfile, mystring, ',');
        dens_vals[i] = std::stod(mystring);
        i++;
    }
    densfile.close();
    y_values = new density(epsilon->get_len(), epsilon, dens_vals);
    delete[] dens_vals;
    
    //setting initial conditions
    x_value = x0;
    dx_value = dx0;
    //NOTE: we never need to use set_ics; setting y_values, x_value, and dx_value happens directly here in the constructor


}

QKESolveMPI::~QKESolveMPI(){
    delete dummy_v_vac;
    for(int i=0; i<epsilon->get_len(); i++){
        delete int_objects[i];
    }
    delete[] int_objects;
    delete just_h;
    delete epsilon;
}
 
void QKESolveMPI::RKCash_Karp(double x, density* y, double dx, double* x_stepped, density* y_5th, density* y_4th)
{
    //int N;
    int N = y->length();
    //to use k1 - need to declare it as an array of N doubles and allocate memory (remember to delete after)
    density* k1 = new density(y);
    density* k2 = new density(y);
    density* k3 = new density(y);
    density* k4 = new density(y);
    density* k5 = new density(y);
    density* k6 = new density(y);
    
    density* z2 = new density(y);
    density* z3 = new density(y);
    density* z4 = new density(y);
    density* z5 = new density(y);
    density* z6 = new density(y); //inputs to get k2-k6
    
    // k1 = dx * f(x, y)
    dx = first_derivative(x, y, k1, dx, x_stepped);
    k1 -> multiply_by(dx);  //k1 = dx * f(x,y)
  
    // k2 = dx * f(x + a2*dx, y + b21*k1)
    z2 -> copy(y);           //z2 = y
    z2 -> add_to(b21, k1);      //z2 = y + b21*k1
    f(x + a2*dx, z2, k2);          //k2 = f(x+a2*dx, z2)
    k2 -> multiply_by(dx);     //dx*f(..)

    //k2->print(8,1);
    // k3 = dx * f(x + a3*dx, y + b31*k1 + b32*k2)
    z3 -> copy(y);           //z3 = y
    z3 -> add_to(b31, k1); //z3 = y + b31*k1
    z3 -> add_to(b32, k2);
    f(x + a3*dx, z3, k3);         // k3 = f(x + a3*dx, z3)
    k3 -> multiply_by(dx);  // k3 = dx*f(x + a3*dx, z3)
 
    // k4 = dx * f(x + a4*dx, y + b41*k1 + b42*k2 +b43*k3)
    z4 -> copy(y);           //z4 = y
    z4 -> add_to(b41, k1);  //z4 = y + b41*k1
    z4 -> add_to(b42, k2); //z4 = y + b41*k1 + b42*k2
    z4 -> add_to(b43, k3); //z4 = y + b41*k1 + b42*k2 + b43*k3
    f(x + a4*dx, z4, k4);         // k4 = f(x + a4*dx, z4)
    k4 -> multiply_by(dx);
        
    // k5 = dx * f(x + a5*dx, y + b51*k1 + b52*k2 + b53*k3 + b54*k4)
    z5 -> copy(y);           //z5 = y
    z5 -> add_to(b51, k1);      //z5 = y + b51*k1
    z5 -> add_to(b52, k2);      //z5 = y + b51*k1 + b52*k2
    z5 -> add_to(b53, k3);      //z5 = y + b51*k1 + b52*k2 + b53*k3
    z5 -> add_to(b54, k4);      //z5 = y + b51*k1 + b52*k2 + b53*k3 + b54*k4    
    f(x + a5*dx, z5, k5);         // k5 = f(x + a5*dx, z5)
    k5 -> multiply_by(dx);
    
    // k6 = dx * f(x + a6*dx, y + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5)
    z6 -> copy(y);           //z6 = y
    z6 -> add_to(b61, k1);      //z6 = y + b61*k1
    z6 -> add_to(b62, k2);      //z6 = y + b61*k1 + b62*k2
    z6 -> add_to(b63, k3);      //z6 = y + b61*k1 + b62*k2 + b63*k3
    z6 -> add_to(b64, k4);      //z6 = y + b61*k1 + b62*k2 + b63*k3 + b64*k4
    z6 -> add_to(b65, k5);      //z6 = y + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5 
    f(x + a6*dx, z6, k6);         // k6 = f(x + a6*dx, z6)
    k6 -> multiply_by(dx);
     
    //y_5th = y + c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5 + c6*k6
    y_5th -> copy(y); //y_5th = y
    y_5th -> add_to(c1, k1);
    y_5th -> add_to(c2, k2);
    y_5th -> add_to(c3, k3);
    y_5th -> add_to(c4, k4);
    y_5th -> add_to(c5, k5);
    y_5th -> add_to(c6, k6);


    // y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3 + cstar4*k4 + cstar5*k5 + cstar6*k6
    y_4th -> copy(y); //y_4th = y           
    y_4th -> add_to(cstar1, k1); //y_4th = y + cstar1*k1
    y_4th -> add_to(cstar2, k2); //y_4th = y + cstar1*k1 + cstar2*k2
    y_4th -> add_to(cstar3, k3); //y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3
    y_4th -> add_to(cstar4, k4); //y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3 + cstar4*k4
    y_4th -> add_to(cstar5, k5); //y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3 + cstar4*k4 + cstar5*k5
    y_4th -> add_to(cstar6, k6); //y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3 + cstar4*k4 + cstar5*k5 + cstar6*k6

    // x_stepped = x + dx
    *x_stepped = x + dx;
    
    delete k1;
    delete k2;
    delete k3;
    delete k4;
    delete k5;
    delete k6;
   
    delete z2;
    delete z3;
    delete z4;
    delete z5;
    delete z6;

    return;
}

bool QKESolveMPI::step_accept(density* y, density* y5, density* y4, double dx, double* dx_new)
{
    //everyone needs to have access to this
    //do i just have one processor do the work and then broadcast the results?
    
    bool accept;
    if (myid == 0){
        int N = y->length();

        int problem = 0;

        double dsm = 0;
        double delta1 = 0;
        double delta0 = 0;

        for (int i = 0; i<N; i++)
        { 
            delta1 = abs(y5 -> get_value(i) - y4 -> get_value(i));
            delta0 = eps*(abs(y -> get_value(i)) + abs(y5 -> get_value(i) - y -> get_value(i))) + TINY;

            if (delta1/delta0 > dsm)
            { 
                dsm = delta1/delta0;
                problem = i;

             }
         }

        if (dsm == 0)
        {
            *dx_new = 5 * dx;
            //cout<< "TRUE (dsm == 0) dx_new = " << *dx_new << endl;
            accept = true;
        } 
        else if (dsm < 1){
            *dx_new = Safety * dx * pow(dsm, -0.2);
            *dx_new = std::min(5.0 * dx, *dx_new); 
            //cout<< "TRUE (dsm < 1) dx_new = " << *dx_new << endl;
            accept = true;
        }
        else{
            *dx_new = Safety * dx * pow(dsm, -0.25);
            //cout<< "FALSE dx_new = " << *dx_new << ", dsm = " << dsm << "; dx = " << dx << endl;
            //cout<< "    i= " << problem << "; y5 = " << y5->get_value(problem) << "; y4 = " << y4->get_value(problem) << endl;

            accept = false;
        }
    }
    MPI_Bcast(&accept, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    return accept;
    
    
}

bool QKESolveMPI::RKCK_step(double x, density* y, double dx, double* x_next, density* y_next, double* dx_next)
{
    //everyone does everything except print out the error message
    double dx_try = dx;
    int N = y->length();
    //dep* y5(y);
    //dep* y4(y);
    density* y5 = new density(y); //???
    density* y4 = new density(y);
    
    bool accept = false;
    for (int i = 0; i<10; i++){
        RKCash_Karp(x, y, dx_try, x_next, y5, y4);
        
        if (step_accept(y, y5, y4, dx_try, dx_next)){
            y_next -> copy(y5);
            accept = true;
            break;
        } 
        else{
           dx_try = *dx_next; 
        }

    }
    if(myid == 0){
        if (!accept)
        {
            cout << "ERROR:  10 iterations without acceptable step" << endl;
            cout << "x = " << x << "; dx = " << dx_try << endl;
        }
        
    }

    delete y5;
    delete y4;
    return accept;
    
    
}

bool QKESolveMPI::ODEOneRun(double x0, density* y0, double dx0, int N_step, int dN, double x_final, double* x, density* y, double* dx, const std::string& file_name, bool verbose) 
{
    //everyone does everything except write the results to the file
    // Set x, y, dx to initial values
    int N = y -> length();
    *x = x0;
    y -> copy(y0);
    *dx = dx0;

    // Declare for RKCK_step
    double* x_next = new double; 
    density* y_next = new density(y);
    double* dx_next = new double; 

    bool no_error = true;
    bool done = false;
    
    if(myid == 0){
        ofstream file(file_name);
        auto start = high_resolution_clock::now();
        
        if (verbose){
            cout << "*******************" << endl;
            cout << "Running ODE Solver.  Initial Conditions:" << endl;
            print_state();
            cout << "Output printed to " << file_name << endl;
        }

        print_csv(file, *x, *dx, y);

        for (int i = 0; i < N_step && no_error && !done; i++){
            for (int j = 0; j < dN; j++){
               // cout << "ith step: " << i << ", jth step: " << j << endl;
                if (*x + *dx > x_final){
                    *dx = x_final - *x;
                }
                
                if (RKCK_step(*x, y, *dx, x_next, y_next, dx_next)){
                    // Update x, y, dx with the results from the RKCK step
                   // cout << "Before update... x =  " << *x << "and dx = " << *dx << endl;
                    *x = *x_next;
                    y->copy(y_next);
                    *dx = *dx_next;
                   // cout << "After update... x = " << *x << "and dx = " << *dx << endl;

                } 
                else{
                    //delete x_next;
                    //delete y_next;
                    //delete dx_next;
                    //file.close();
                    no_error = false;
                    break;
                    //return false;
                }

                if (*x == x_final){
                    cout << "Reached x_final" << endl;
                    print_csv(file, *x, *dx, y);
                    //delete x_next;
                    //delete y_next;
                    //delete dx_next;
                    //file.close();
                    done = true;
                    break;
                    //return true;
                }
            }

            print_csv(file, *x, *dx, y_next);

        }

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);

        if (verbose){
            print_state();
            cout << endl << "Time elapsed: "
             << duration.count()/1000. << " seconds" << endl;
        }
        file.close();
        
    }
    
    else{
        for (int i = 0; i < N_step && no_error && !done; i++){
            for (int j = 0; j < dN; j++){
                if (*x + *dx > x_final){
                    *dx = x_final - *x;
                }
                
                if (RKCK_step(*x, y, *dx, x_next, y_next, dx_next)){
                    *x = *x_next;
                    y->copy(y_next);
                    *dx = *dx_next;
                } 
                
                else{
                    no_error = false;
                    break;
                }

                if (*x == x_final){
                    done = true;
                    break;
                }
            }
        
        }
    }

    delete x_next;
    delete y_next;
    delete dx_next; 
    
    return true;
}

//no need to overwrite print_csv or print_state as these functions should only ever be called by main


bool QKESolveMPI::run(int N_step, int dN, double x_final, const std::string& file_name, bool verbose)
{
    //we want everyone to run this
    return ODEOneRun(x_value, y_values, dx_value, N_step, dN, x_final, &x_value, y_values, &dx_value, file_name, verbose);

}

void QKESolveMPI::f(double t, density* d1, density* d2)
{
    d2->zeros();
    double* d2_vals = new double[d1->length()]();
    double myans=0;
    int sender, tag;
    double* dummy_int = new double[4];
    
    MPI_Status status;
        
    
    if(myid == 0){
        three_vector_for_QKE* dummy_v_dens = new three_vector_for_QKE;
        three_vector_for_QKE* dummy_v_therm = new three_vector_for_QKE;

        three_vector* V_nu = new three_vector;
        three_vector* V_nubar = new three_vector;
        three_vector* p = new three_vector;
        three_vector* vcrossp = new three_vector;

        dummy_v_dens->v_density(epsilon, d1);
        dummy_v_therm->v_thermal(epsilon, d1);
        double Tcm = d1->get_Tcm();
        double en = 0.;
        for (int i=1; i< epsilon->get_len(); i++){
            en = epsilon->get_value(i) * Tcm;
            
            V_nu->copy(dummy_v_dens);
            V_nu->add_to(1./en, dummy_v_vac);
            V_nu->add_to(en, dummy_v_therm);

            d1->p_vector(i, true, p);

            vcrossp->set_cross_product(V_nu, p);
            d2_vals[4*i+1] = vcrossp->get_value(0);
            d2_vals[4*i+2] = vcrossp->get_value(1);
            d2_vals[4*i+3] = vcrossp->get_value(2);

            V_nubar->copy(dummy_v_dens);
            V_nubar->add_to(-1./en, dummy_v_vac);
            V_nubar->add_to(-en, dummy_v_therm);

            d1->p_vector(i, false, p);
            vcrossp->set_cross_product(V_nubar, p);
            d2_vals[4*(epsilon->get_len())+4*i+1] = vcrossp->get_value(0);
            d2_vals[4*(epsilon->get_len())+4*i+2] = vcrossp->get_value(1);
            d2_vals[4*(epsilon->get_len())+4*i+3] = vcrossp->get_value(2);
            
        }
        for(int i=0; i<2*epsilon->get_len(); i++){
            MPI_Recv(dummy_int, 4, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            sender = status.MPI_SOURCE;
            tag = status.MPI_TAG;
            for(int j=0; j<4; j++){
                d2_vals[4*tag+j] += dummy_int[j];
            }
        }
        
        
        delete dummy_v_dens;
        delete dummy_v_therm;
        delete V_nu;
        delete V_nubar;
        delete vcrossp;
        delete p;
        
        
    }
    
    else{
        //OTHER PROCESSORS FIND INTEGRALS AND SEND BACK TO MAIN
        for(int i=myid-1; i<epsilon->get_len(); i+=numprocs-1){
            
            //antineutrino 
            int_objects[i]->whole_integral(d1, false, dummy_int);
            MPI_Send(dummy_int, 4, MPI_DOUBLE, 0, epsilon->get_len()+i, MPI_COMM_WORLD);
            
            //neutrino
            int_objects[i]->whole_integral(d1, true, dummy_int);
            MPI_Send(dummy_int, 4, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
        }
        
    } 
    //MAIN BROADCASTS OUT D2 AS A VALUES ARRAY, EVERYONE RECIEVES AND CONVERTS TO DENSITY OBJECT
    MPI_Bcast(d2_vals, d1->length(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for(int i=0; i<d1->length(); i++){
        d2->set_value(i, d2_vals[i]);
    }
    
    delete[] dummy_int;
    delete[] d2_vals;

}

double QKESolveMPI::first_derivative(double t, density* d1, density* d2, double dx, double* x_next){
    d2->zeros();
    double* d2_vals = new double[d1->length()]();
    double myans = 0;
    double new_dx = 0;
    int sender, tag;
    double* dummy_int = new double[4];
    
    MPI_Status status;
        
    
    if(myid == 0){
        double* next_dx = new double;
        density* y_next = new density(d1);
        just_h->RKCK_step(t, d1, dx, x_next, y_next, next_dx);
        new_dx = *x_next-t;
        delete next_dx;
        delete y_next;
        
        three_vector_for_QKE* dummy_v_dens = new three_vector_for_QKE;
        three_vector_for_QKE* dummy_v_therm = new three_vector_for_QKE;

        three_vector* V_nu = new three_vector;
        three_vector* V_nubar = new three_vector;
        three_vector* p = new three_vector;
        three_vector* vcrossp = new three_vector;

        dummy_v_dens->v_density(epsilon, d1);
        dummy_v_therm->v_thermal(epsilon, d1);
        double Tcm = d1->get_Tcm();
        double en = 0.;
        for (int i=1; i< epsilon->get_len(); i++){
            en = epsilon->get_value(i) * Tcm;
            
            V_nu->copy(dummy_v_dens);
            V_nu->add_to(1./en, dummy_v_vac);
            V_nu->add_to(en, dummy_v_therm);

            d1->p_vector(i, true, p);

            vcrossp->set_cross_product(V_nu, p);
            d2_vals[4*i+1] = vcrossp->get_value(0);
            d2_vals[4*i+2] = vcrossp->get_value(1);
            d2_vals[4*i+3] = vcrossp->get_value(2);

            V_nubar->copy(dummy_v_dens);
            V_nubar->add_to(-1./en, dummy_v_vac);
            V_nubar->add_to(-en, dummy_v_therm);

            d1->p_vector(i, false, p);
            vcrossp->set_cross_product(V_nubar, p);
            d2_vals[4*(epsilon->get_len())+4*i+1] = vcrossp->get_value(0);
            d2_vals[4*(epsilon->get_len())+4*i+2] = vcrossp->get_value(1);
            d2_vals[4*(epsilon->get_len())+4*i+3] = vcrossp->get_value(2);
            
        }
        
        for(int i=0; i<2*epsilon->get_len(); i++){
            MPI_Recv(dummy_int, 4, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            sender = status.MPI_SOURCE;
            tag = status.MPI_TAG;
            for(int j=0; j<4; j++){
                d2_vals[4*tag+j] += dummy_int[j];
            }
        }
        
        delete dummy_v_dens;
        delete dummy_v_therm;
        delete V_nu;
        delete V_nubar;
        delete vcrossp;
        delete p;
        
        
    }
    
    else{
        //OTHER PROCESSORS FIND INTEGRALS AND SEND BACK TO MAIN
        
        
        for(int i=myid-1; i<epsilon->get_len(); i+=numprocs-1){
            
            //antineutrino 
            int_objects[i]->whole_integral(d1, false, dummy_int);
            MPI_Send(dummy_int, 4, MPI_DOUBLE, 0, epsilon->get_len()+i, MPI_COMM_WORLD);
            
            //neutrino
            int_objects[i]->whole_integral(d1, true, dummy_int);
            MPI_Send(dummy_int, 4, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
        }
        
    } 
    //MAIN BROADCASTS OUT D2 AS A VALUES ARRAY, EVERYONE RECIEVES AND CONVERTS TO DENSITY OBJECT
    MPI_Bcast(d2_vals, d1->length(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for(int i=0; i<d1->length(); i++){
        d2->set_value(i, d2_vals[i]);
    }
    MPI_Bcast(&new_dx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    delete[] dummy_int;
    return new_dx;
}
