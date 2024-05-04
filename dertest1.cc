// import statements
#include <iostream>
#include <cmath>
#include <algorithm>
#include "constants.hh" // so now I can just write the names of variables in constants.hh
#include "CashKarp_vals.hh"
#include "arrays.hh"



using std::cout;
using std::endl;
/*
using std::pow
using std::cos
using std::sin
do i have to do this if im just saying cos() later???????????????????????????????????????????????????
*/

//void f(double, double*, int, double*);
//void RKCash_Karp(double, double*, double, int, double*, double*, double*);
//void scalar_times_vector(double, double*, int);
//void copy_vector(double*, double*, int);
//void add_vector(double, double*, double*, int);


void f(double, dep_vars*, dep_vars*);
void RKCash_Karp(double, dep_vars*, double, double*, dep_vars*, dep_vars*);
bool step_accept(dep_vars*, dep_vars*, dep_vars*, double, double*);
bool RKCK_step(double, dep_vars*, double, double*, dep_vars*, double*);
//void dep_vars::multiply_by(double);
//void dep_vars::copy(dep_vars*);
//void dep_vars::add_to(double, dep_vars*);








int main(int argc, char** argv) 
{     
    int N;
    if (argc == 5) 
    {
        N = 4;
        //cout << "Using der2 " << endl;
        
    } else { // NOTE: for some reason you have to put the } else { on the same line
        //cout << "Using der1 " << endl;
        N = 3;
    }
    
/******************************************
SIDE NOTE:

argv is an "array of C-style strings (char*) representing the actual arguments passed to the program"
- argv[0] is the name of the program
- argv[1], argv[2], etc., are the following arguments you put in

argc is the count of how many things (strings) are in argv (including name of program)
- so here this is saying if argc > 5, or if there are 5 arrays in argv??? use N=4
- And if argc is not larger than 5, use N=3
- changing from argc changes how many function evaluations and derivative calculations are performed
- but its about like how many things you put in the terminal?

The logic here is that in the terminal when i run the program "./testder" if i give it 1-3 arguments (so like ./testder is 1 argument, ./testder 0.75 is 2 arguments, etc.) it will run der1 and N=3. However, if i give it 5 arguments so the name of the program and then 4 input values of y, it will run der2 the one with 4 dependent variables. The number of arguments given is argc so thats why the if else statement works this way. 
*******************************************/
    
    //initialize variables
    double x = 0;
    double x_stepped = 0;
    double dx = 0.1;
    
   // double* y = new double[N]();
   // double* y_5th = new double[N]();
   // double* y_4th = new double[N]();
   // double *der = new double[N];
        //these become:
    dep_vars* y = new dep_vars(N);
    dep_vars* y_5th = new dep_vars(N);
    dep_vars* y_4th = new dep_vars(N);
    dep_vars* der = new dep_vars(N);
    

    if (argc == 5)
    {       
       for (int i = 0; i<N; i++)
       {
           y -> set_value(i, atof(argv[i+1]));
       }
     }
    
    f(x, y, der);
    
    RKCash_Karp(x, y, dx, &x_stepped, y_5th, y_4th); // the & delivers the memory address , a pointer
    
    //print y5th and y4th
    y_5th -> print_all();
    y_4th -> print_all();

    
    // delete them
    delete y;
    delete der;
    delete y_5th;
    delete y_4th;
    
    return 0; 
}














void f(double x, dep_vars* y, dep_vars* z)
{
    int N = y->length();

    if (N==3)
    {
        z->set_value(0, 1.);
        z->set_value(1, x);
        z->set_value(2, 1/(1+x));
    }
    else if (N==4)
    {
        z->set_value(0, y->get_value(1));
        z->set_value(1, y->get_value(2));
        z->set_value(2, y->get_value(3));
        z->set_value(3, pow(_PI_, 4) * y->get_value(0));
    }
    return;
}









/***************************
double* x_stepped is a pointer to the memory spot of a single double precision number
double* y_5th and 4th are pointers to arrays of doubles each with N elements
***************************/
void RKCash_Karp(double x, dep_vars* y, double dx, double* x_stepped, dep_vars* y_5th, dep_vars* y_4th)
{
    //int N;
    int N = y->length();
    //to use k1 - need to delcaire it as an array of N doubles and allocate memory (remember to delete after)
    dep_vars* k1 = new dep_vars(N);
    dep_vars* k2 = new dep_vars(N);
    dep_vars* k3 = new dep_vars(N);
    dep_vars* k4 = new dep_vars(N);
    dep_vars* k5 = new dep_vars(N);
    dep_vars* k6 = new dep_vars(N);
    
    dep_vars* z2 = new dep_vars(N);
    dep_vars* z3 = new dep_vars(N);
    dep_vars* z4 = new dep_vars(N);
    dep_vars* z5 = new dep_vars(N);
    dep_vars* z6 = new dep_vars(N); //inputs to get k2-k6
    
//k1
    // k1 = dx * f(x, y)
    f(x, y, k1);
    k1 -> multiply_by(dx);  //k1 = dx * f(x,y)
    //cout << "k1: ";
   // for(int i = 0; i < N; i++) cout << k1[i] << " ";
    //cout << endl;
    
//k2
    
    // k2 = dx * f(x + a2*dx, y + b21*k1)
    z2 -> copy(y);           //z2 = y
    z2 -> add_to(b21, k1);      //z2 = y + b21*k1
    f(x + a2*dx, z2, k2);          //k2 = f(x+a2*dx, z2)
    k2 -> multiply_by(dx);     //dx*f(..)
    //cout << "k2: ";
    //for(int i = 0; i < N; i++) cout << k2[i] << " ";
    //cout << endl;
    
    
    
//k3
    
    // k3 = dx * f(x + a3*dx, y + b31*k1 + b32*k2)
    z3 -> copy(y);           //z3 = y
    //cout << "z3 = y ----> ";
    //for(int i = 0; i < N; i++) cout << z3[i] << " ";
    //cout << endl;
    //add_vector(b31, k1, z3, N);      //z3 = y + b31*k1
    z3 -> add_to(b31, k1); //z3 = y + b31*k1
   // cout << "z3 = y + b31*k1 ----> ";
    //for(int i = 0; i < N; i++) cout << z3[i] << " ";
    //cout << endl; 
    //add_vector(b32, k2, z3, N);      //z3 = y + b31*k1 + b32*k2
    z3 -> add_to(b32, k2);
    //cout << "z3 = y + b31*k1 + b32*k2 ---->";
    //for(int i = 0; i < N; i++) cout << z3[i] << " ";
    //cout << endl;
    f(x + a3*dx, z3, k3);         // k3 = f(x + a3*dx, z3)
    //cout << "k3: ";
    //for(int i = 0; i < N; i++) cout << k3[i] << " ";
    //cout << endl;
    k3 -> multiply_by(dx);  // k3 = dx*f(x + a3*dx, z3)
   // cout << "k3: ";
    //for(int i = 0; i < N; i++) cout << k3[i] << " ";
    //cout << endl;
    
    
//k4
        
    // k4 = dx * f(x + a4*dx, y + b41*k1 + b42*k2 +b43*k3)
    z4 -> copy(y);           //z4 = y
    z4 -> add_to(b41, k1);  //z4 = y + b41*k1
    z4 -> add_to(b42, k2); //z4 = y + b41*k1 + b42*k2
    z4 -> add_to(b43, k3); //z4 = y + b41*k1 + b42*k2 + b43*k3
    f(x + a4*dx, z4, k4);         // k4 = f(x + a4*dx, z4)
    k4 -> multiply_by(dx);
    //cout << "k4: ";
    //for(int i = 0; i < N; i++) cout << k4[i] << " ";
    //cout << endl;
        
    // k5 = dx * f(x + a5*dx, y + b51*k1 + b52*k2 + b53*k3 + b54*k4)
    z5 -> copy(y);           //z5 = y
    z5 -> add_to(b51, k1);      //z5 = y + b51*k1
    z5 -> add_to(b52, k2);      //z5 = y + b51*k1 + b52*k2
    z5 -> add_to(b53, k3);      //z5 = y + b51*k1 + b52*k2 + b53*k3
    z5 -> add_to(b54, k4);      //z5 = y + b51*k1 + b52*k2 + b53*k3 + b54*k4    
    f(x + a5*dx, z5, k5);         // k5 = f(x + a5*dx, z5)
    k5 -> multiply_by(dx);
    //cout << "k5: ";
    //for(int i = 0; i < N; i++) cout << k5[i] << " ";
    //cout << endl;
    
    // k6 = dx * f(x + a6*dx, y + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5)
    z6 -> copy(y);           //z6 = y
    z6 -> add_to(b61, k1);      //z6 = y + b61*k1
    z6 -> add_to(b62, k2);      //z6 = y + b61*k1 + b62*k2
    z6 -> add_to(b63, k3);      //z6 = y + b61*k1 + b62*k2 + b63*k3
    z6 -> add_to(b64, k4);      //z6 = y + b61*k1 + b62*k2 + b63*k3 + b64*k4
    z6 -> add_to(b65, k5);      //z6 = y + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5 
    f(x + a6*dx, z6, k6);         // k6 = f(x + a6*dx, z6)
    k6 -> multiply_by(dx);
    
   // cout << "k6: ";
    //for(int i = 0; i < N; i++) cout << k6[i] << " ";
    //cout << endl;
    
        
    //y_5th = y + c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5 + c6*k6
    y_5th -> copy(y); //y_5th = y
   // add_vector(c1, k1, y_5th, N);    //y_5th = y + c1*k1
    y_5th -> add_to(c1, k1);
   // add_vector(c2, k2, y_5th, N);    //y_5th = y + c1*k1 + c2*k2
    y_5th -> add_to(c2, k2);
   // add_vector(c3, k3, y_5th, N);    //y_5th = y + c1*k1 + c2*k2 + c3*k3
    y_5th -> add_to(c3, k3);
   // add_vector(c4, k4, y_5th, N);    //y_5th = y + c1*k1 + c2*k2 + c3*k3 + c4*k4
    y_5th -> add_to(c4, k4);
   // add_vector(c5, k5, y_5th, N);    //y_5th = y + c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5
    y_5th -> add_to(c5, k5);
   // add_vector(c6, k6, y_5th, N);    //y_5th = y + c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5 + c6*k6
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


















double eps = 1e-8;
double TINY = 1e-40;
double Safety = 0.9;
    
bool step_accept(dep_vars* y, dep_vars* y5, dep_vars* y4, double dx, double* dx_new)
{
    int N = y->length();
    
    double dsm = 0;
    double delta1;
    double delta0;
    
    for (int i = 0; i<N; i++)
    { 
        delta1 = abs(y5 -> get_value(i) - y4 -> get_value(i));
        delta0 = eps*(abs(y -> get_value(i)) + abs(y5 -> get_value(i) - y -> get_value(i)));
        
        if (delta1/delta0 > dsm)
        { 
            dsm = delta1/delta0;
         }
     }
     
    if (dsm == 0)
    {
        *dx_new = 5 * dx;
        return true;
    } 
    else if (dsm <1){
        *dx_new = Safety * dx * pow(dsm, -0.25);
        *dx_new = std::min(5.0 * dx, *dx_new); 
        return true;
    }
    else{
        *dx_new = Safety * dx * pow(dsm, -0.25);
        return false;
    }
    
    
}














bool RKCK_step(double x, dep_vars* y, double dx, double* x_next, dep_vars* y_next, double* dx_next)
{
    double dx_try = dx;
    int N = y->length();
    dep_vars* y5 = new dep_vars(N); //???
    dep_vars* y4 = new dep_vars(N);
    
    for (int i = 0; i<10; i++)
    {
        RKCash_Karp(x, y, dx_try, x_next, y5, y4);
        if (step_accept(y, y5, y4, dx_try, dx_next))
        {
            y_next -> copy(y5);
            return true;
            break;
        } 
        else {
           dx_try = *dx_next; 
        }
        
    }
    
    cout << "ERROR:  10 iterations without acceptable step" << endl;
    
    delete y5;
    delete y4;
    
    return false;
    
    
}












    