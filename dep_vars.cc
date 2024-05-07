#include <iostream>
#include "arrays.hh"

using std::cout;
using std::endl;
using std::ostream;

dep_vars::dep_vars(int size)
{
    N = size;
    values = new double[N]();
}
    
dep_vars::dep_vars(double* copy_me, int size)
{
    N = size;
    values = new double[N]();
    for (int i = 0; i < N; i++)
        values[i] = copy_me[i];
}
    
dep_vars::dep_vars(dep_vars* copy_me)
{
    N = copy_me->length();
    values = new double[N]();
    for (int i = 0; i < N; i++)
        values[i] = copy_me->values[i];
}
    
dep_vars::~dep_vars()
{delete[] values;}

int dep_vars::length()
{return N;}

double dep_vars::get_value(int i)
{return values[i];}

void dep_vars::set_value(int i, double v)
{values[i] = v;}
    
void dep_vars::print_all()
{
    for (int i = 0; i < N; i++)
        cout << values[i] << endl;
}
    
void dep_vars::print(int N_top = 3, int N_bot = 1)
{
    if (N <= N_top + N_bot)
        print_all();
    else if (N_top < 0 || N_bot < 0)
        print_all();
    else
    {
        for (int i = 0; i < N_top; i++)
            cout << values[i] << endl;
        cout << "..." << endl;
        for (int i = 0; i < N_bot; i++)
            cout << values[N - N_bot + i] << endl;
    }
}

/********************************
/ We need to write the tested methods:
/ void dep_vars::copy(dep_vars*)
/ void dep_vars::multiply_by(double)
/ void dep_vars::add_to(double, dep_vars*)

/================================
/ These relate to: 
/   void copy_vector(double* z, double* y, int N)
/   void scalar_times_vector(double c, double* y, int N)
/   void add_vector(double c, double* z, double* y, int N)
/
/ in that (double* y, int N) are information stored as a part of the dep_vars class
/ as a result, instead of writing y[i] we will write values[i] and we can use N without it being an input
/ double* z is replaced by dep_vars* z, another instance of the dep_vars class.
/ instead of writing z[i], we will use z->get_value(i)
/ 
/ This will be a cleaner version of the old functions, and protects the values within the functions
***********************************/

void dep_vars::multiply_by(double scalar)
{
    for (int i = 0; i < N; ++i) 
    {
        values[i] *= scalar;
    }
}






void dep_vars::copy(dep_vars* z)
{
    
    for (int i = 0; i < N; ++i) 
    {
        values[i] = z -> get_value(i);
    }
}



void dep_vars::add_to(double c, dep_vars* z)
{

    for (int i = 0; i < N; ++i) {
        //values[i] += c*z[i];
        values[i] += c * z -> get_value(i);
    }

}


        
