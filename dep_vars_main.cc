#include <iostream>
#include <cmath>
#include "arrays.hh"
#include "constants.hh" // so now I can just write the names of variables in constants.hh


using std::cout;
using std::endl;

/*
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
*/

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

int main()
{
    double x = 0.5;
    double y[] = {1.1, 2.2, 3.3};
    dep_vars* input = new dep_vars(y, 3);
    dep_vars* output = new dep_vars(3);
    
    f(x, input, output);
    
    output->print_all();

    three_vector* tv = new three_vector(1, 2, 1);
    cout << tv->magnitude_squared() << endl;

    three_vector* tv2 = new three_vector(tv);
    cout << tv2->dot_with(tv) << endl;
    
    delete tv;
    delete tv2;
    
    delete input;
    delete output;
    return 0;
}