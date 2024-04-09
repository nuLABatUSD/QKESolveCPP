// import statements
#include <iostream>
using std::cout;
using std::endl;

// have to put the signature bc i defined the function later
void f(double, double*, int, double*);


//main funciton
int main(int argc, char** argv) { 
    double x = 0;
    if (argc == 2){
        x = atof(argv[1]);}
                 
    int N = 3;
    double* y = new double[N]();
    double *der = new double[N];
    
    
    // call the function
    f(x, y, N, der);

    // print out the results for der
    cout << "Results of der:" << endl;
    for (int i = 0; i < N; ++i) {
        cout << "der[" << i << "] = " << der[i] << endl;}
    
    // delete them
    delete[] y;
    delete[] der;
    
    return 0; }

//functions
void f(double x, double* y, int N, double* der){
    der[0] = 1;
    der[1] = x;
    der[2] = 1/(x+1);}

