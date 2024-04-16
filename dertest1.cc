// import statements
#include <iostream>
#include <cmath>
#include "constants.hh" // so now I can just write the names of variables in constants.hh
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
    double* y = new double[N]();
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

    // print out the results for der
    cout << "Results of der and y values:" << endl;
    for (int i = 0; i < N; ++i) 
    {
        cout << "der[" << i << "] = " << der[i] << endl;
        cout << "y[" << i << "] = " << y[i] << endl;
    }
    
    // delete them
    delete[] y;
    delete[] der;
    
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
// WHY DO I USE ELSE IF INSTEAD OF ELSE HERE??????????????????????????????????????????????????
    } else if (N == 4) {
        der[0] = y[1];
        der[1] = y[2];
        der[2] = y[3];
        der[3] = pow(_PI_, 4) * y[0];
    } // really its std::pow?
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
   