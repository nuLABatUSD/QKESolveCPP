// import statements
#include <iostream>
#include <cmath>
#include "constants.hh" // so now I can just write the names of variables in constants.hh
#include "CashKarp_vals.hh"

using std::cout;
using std::endl;
/*
using std::pow
using std::cos
using std::sin
do i have to do this if im just saying cos() later???????????????????????????????????????????????????
*/

// have to put the signature bc i defined the function later
void f(double, double*, int, double*);
void RKCash_Karp(double, double*, double, int, double*, double*, double*);
void scalar_times_vector(double, double*, int);
void copy_vector(double*, double*, int);
void add_vector(double, double*, double*, int);








//main funciton
int main(int argc, char** argv) 
{     
    int N;
    if (argc == 5) 
    {
        N = 4;
        cout << "Using der2 " << endl;
        
    } else { // NOTE: for some reason you have to put the } else { on the same line
        cout << "Using der1 " << endl;
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
    double* y = new double[N]();
    double* y_5th = new double[N]();
    double* y_4th = new double[N]();
    double *der = new double[N];

    if (argc == 5)
    {       
       for (int i = 0; i<N; i++)
       {
            y[i] = atof(argv[i + 1]); //this says: y[0] = argv[0], y[1] = argv[1] and so on. this means that whatever i input in the terminal when i call the function for my arguments will be filled in as the y values! In doing so, the value given for argument [i] is converted to a double (thats the atof)?
       }
     }
    
    // call the function
    f(x, y, N, der);
    
    RKCash_Karp(x, y, dx, N, &x_stepped, y_5th, y_4th); // the & delivers the memory address , a pointer
    
    //print x_stepped
    cout << "x_stepped = " << x_stepped << endl;
    
    //print y5th and y4th
     for (int i = 0; i < N; ++i) 
    {
        cout << "y_5th = " << y_5th[i] << endl;
    }
    
         for (int i = 0; i < N; ++i) 
    {
        cout << "y_4th = " << y_4th[i] << endl;
    }
    
    // delete them
    delete[] y;
    delete[] der;
    delete[] y_5th;
    delete[] y_4th;
    
    return 0; 
}











//functions
void f(double x, double* y, int N, double* der) 
{
    if (N == 3)
    {
        der[0] = 1;
        der[1] = x;
        der[2] = 1/(x+1);

    } else if (N == 4) {
        der[0] = y[1];
        der[1] = y[2];
        der[2] = y[3];
        der[3] = pow(_PI_, 4) * y[0]; // really its std::pow?
    } 
} 










/***************************
double* x_stepped is a pointer to the memory spot of a single double precision number
double* y_5th and 4th are pointers to arrays of doubles each with N elements
***************************/
void RKCash_Karp(double x, double* y, double dx, int N, double* x_stepped, double* y_5th, double* y_4th)
{
    
    //to use k1 - need to delcaire it as an array of N doubles and allocate memory (remember to delete after)
    double* k1 = new double[N]();
    double* k2 = new double[N]();
    double* k3 = new double[N]();
    double* k4 = new double[N]();
    double* k5 = new double[N]();
    double* k6 = new double[N]();
    
    double* z2 = new double[N]();
    double* z3 = new double[N]();
    double* z4 = new double[N]();
    double* z5 = new double[N]();
    double* z6 = new double[N](); //inputs to get k2-k6
    
    // k1 = dx * f(x, y)
    f(x, y, N, k1);
    scalar_times_vector(dx, k1, N);  //k1 = dx * f(x,y)
    
    cout << "k1: ";
    for(int i = 0; i < N; i++) cout << k1[i] << " ";
    cout << endl;
    
    
    // k2 = dx * f(x + a2*dx, y + b21*k1)
    copy_vector(y, z2, N);           //z2 = y
    add_vector(b21, k1, z2, N);      //z2 = y + b21*k1
    f(x + a2*dx, z2, N, k2);          //k2 = f(x+a2*dx, z2)
    scalar_times_vector(dx, k2, N);     //dx*f(..)
    
    cout << "k2: ";
    for(int i = 0; i < N; i++) cout << k2[i] << " ";
    cout << endl;
    
    
    

    
    // k3 = dx * f(x + a3*dx, y + b31*k1 + b32*k2)
    copy_vector(y, z3, N);           //z3 = y
    
    cout << "z3 = y ----> ";
    for(int i = 0; i < N; i++) cout << z3[i] << " ";
    cout << endl;
    
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR STARTS HERE IN CALCULATING Z3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    add_vector(b31, k1, z3, N);      //z3 = y + b31*k1
    
    cout << "z3 = y + b31*k1 ----> ";
    for(int i = 0; i < N; i++) cout << z3[i] << " ";
    cout << endl; 
    
    add_vector(b32, k2, z3, N);      //z3 = y + b31*k1 + b32*k2
    
    cout << "z3 = y + b31*k1 + b32*k2 ---->";
    for(int i = 0; i < N; i++) cout << z3[i] << " ";
    cout << endl;
    
 
    f(x + a3*dx, z3, N, k3);         // k3 = f(x + a3*dx, z3)
    
    cout << "k3: ";
    for(int i = 0; i < N; i++) cout << k3[i] << " ";
    cout << endl;
    
    scalar_times_vector(dx, k3, N);  // k3 = dx*f(x + a3*dx, z3)
    
    cout << "k3: ";
    for(int i = 0; i < N; i++) cout << k3[i] << " ";
    cout << endl;
    
    
////////////DEBUGGING K3 AHHAHHAHHAAHHHHHHHHHHHHHHH///////////
        
    // k4 = dx * f(x + a4*dx, y + b41*k1 + b42*k2 +b43*k3)
    copy_vector(y, z4, N);           //z4 = y
    add_vector(b41, k1, z4, N);      //z4 = y + b41*k1
    add_vector(b42, k2, z4, N);      //z4 = y + b41*k1 + b42*k2
    add_vector(b43, k3, z4, N);      //z4 = y + b41*k1 + b42*k2 + b43*k3
    f(x + a4*dx, z4, N, k4);         // k4 = f(x + a4*dx, z4)
    scalar_times_vector(dx, k4, N);
    
    cout << "k4: ";
    for(int i = 0; i < N; i++) cout << k4[i] << " ";
    cout << endl;
        
    // k5 = dx * f(x + a5*dx, y + b51*k1 + b52*k2 + b53*k3 + b54*k4)
    copy_vector(y, z5, N);           //z5 = y
    add_vector(b51, k1, z5, N);      //z5 = y + b51*k1
    add_vector(b52, k2, z5, N);      //z5 = y + b51*k1 + b52*k2
    add_vector(b53, k3, z5, N);      //z5 = y + b51*k1 + b52*k2 + b53*k3
    add_vector(b54, k4, z5, N);      //z5 = y + b51*k1 + b52*k2 + b53*k3 + b54*k4
    f(x + a5*dx, z5, N, k5);         // k5 = f(x + a5*dx, z5)
    scalar_times_vector(dx, k5, N);
    
    cout << "k5: ";
    for(int i = 0; i < N; i++) cout << k5[i] << " ";
    cout << endl;
    
    // k6 = dx * f(x + a6*dx, y + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5)
    copy_vector(y, z6, N);           //z6 = y
    add_vector(b61, k1, z6, N);      //z6 = y + b61*k1
    add_vector(b62, k2, z6, N);      //z6 = y + b61*k1 + b62*k2
    add_vector(b63, k3, z6, N);      //z6 = y + b61*k1 + b62*k2 + b63*k3
    add_vector(b64, k4, z6, N);      //z6 = y + b61*k1 + b62*k2 + b63*k3 + b64*k4
    add_vector(b65, k5, z6, N);      //z6 = y + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5
    f(x + a6*dx, z6, N, k6);         // k6 = f(x + a6*dx, z6)
    scalar_times_vector(dx, k6, N);
    
    cout << "k6: ";
    for(int i = 0; i < N; i++) cout << k6[i] << " ";
    cout << endl;
    
        
    //y_5th = y + c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5 + c6*k6
    copy_vector(y, y_5th, N); //y_5th = y
    add_vector(c1, k1, y_5th, N);    //y_5th = y + c1*k1
    add_vector(c2, k2, y_5th, N);    //y_5th = y + c1*k1 + c2*k2
    add_vector(c3, k3, y_5th, N);    //y_5th = y + c1*k1 + c2*k2 + c3*k3
    add_vector(c4, k4, y_5th, N);    //y_5th = y + c1*k1 + c2*k2 + c3*k3 + c4*k4
    add_vector(c5, k5, y_5th, N);    //y_5th = y + c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5
    add_vector(c6, k6, y_5th, N);    //y_5th = y + c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5 + c6*k6
    
        
    // y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3 + cstar4*k4 + cstar5*k5 + cstar6*k6
    copy_vector(y, y_4th, N); //y_4th = y
    add_vector(cstar1, k1, y_4th, N); //y_4th = y + cstar1*k1
    add_vector(cstar2, k2, y_4th, N); //y_4th = y + cstar1*k1 + cstar2*k2
    add_vector(cstar3, k3, y_4th, N); //y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3
    add_vector(cstar4, k4, y_4th, N); //y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3 + cstar4*k4
    add_vector(cstar5, k5, y_4th, N); //y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3 + cstar4*k4 + cstar5*k5
    add_vector(cstar6, k6, y_4th, N); //y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3 + cstar4*k4 + cstar5*k5 + cstar6*k6
        
    // x_stepped = x + dx
    *x_stepped = x + dx;
  
        
    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
    delete[] k5;
    delete[] k6;
   
    delete[] z2;
    delete[] z3;
    delete[] z4;
    delete[] z5;
    delete[] z6;
       
    return;
}












// going to do a lot of dx * f() but we cant just straight up do that in python so making a function
// want this to take every value in y array and multiply it by c
void scalar_times_vector(double scalar, double* vector, int N)
{
    for (int i = 0; i < N; ++i) 
    {
        vector[i] *= scalar;
    }
}







void copy_vector(double* z, double* y, int N)
{
    
    for (int i = 0; i < N; ++i) 
    {
        y[i] = z[i];
    }
}




void add_vector(double c, double* z, double* y, int N)
{
    double* zn = new double[N];

    copy_vector(z, zn, N); //added this part and it works!

    scalar_times_vector(c, zn, N);

    for (int i = 0; i < N; ++i) {
        y[i] += zn[i];
    }

    delete[] zn;
}





/*
dont actually need the y functions
    //initialize y for N=3
    y[0] = x + 1;
    y[1] = 0.5 * pow(x, 2) + 1;
    y[2] = log(x + 1) + 1; 
    
    //initialize y for N=4
    y[0] = cos(_PI_ * x);
    y[1] = -_PI_ * sin(_PI_ * x);
    y[2] = -pow(_PI_, 2) * cos(_PI_ * x);
    y[3] = pow(_PI_, 3) * sin(_PI_ * x);
*/
   