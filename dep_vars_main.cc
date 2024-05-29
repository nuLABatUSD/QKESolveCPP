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
   
   
   linspace* et = new linspace(0.,20, 201);
   three_vector* v = new three_vector();
   

   
   density* den = new density(et, .01, -0.01);
   cout << den->num_bins() << endl;
   den->p0_p(1, true, v);
   
   cout << "now the stuff i care about: " << endl;
   v->print_all();
   den->print_all();
   

   
   delete et;
   delete den;
   delete v;
   return 0;
}