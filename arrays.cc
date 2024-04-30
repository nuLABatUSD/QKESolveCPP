#include "arrays.hh"
#include <cmath>

three_vector::three_vector():dep_vars(3)
{;}

three_vector::three_vector(double x, double y, double z):dep_vars(3)
{
    values[0] = x;
    values[1] = y;
    values[2] = z;
}

three_vector::three_vector(double* copy_me):dep_vars(copy_me, 3)
{;}

three_vector::three_vector(three_vector* copy_me):dep_vars(copy_me)
{;}

double three_vector::dot_with(three_vector* B)
{
    double dot = 0;
    for(int i = 0; i < 3; i++)
        dot += values[i] * B->get_value(i);
    return dot;
}

double three_vector::magnitude_squared()
{
    return dot_with(this);
}

double three_vector::magnitude()
{
    double sum = 0;
    for(int i =0; i < 3; i++)
        sum += pow(this->get_value(i),2);
    return sqrt(sum);
}

void three_vector::set_cross_product(three_vector* A, three_vector* B)
{
    values[0] = A->get_value(1) * B->get_value(2) - A->get_value(2) * B->get_value(1);
    values[1] = A->get_value(2) * B->get_value(0) - A->get_value(0) * B->get_value(2);
    values[2] = A->get_value(0) * B->get_value(1) - A->get_value(1) * B->get_value(0);
}


linspace::linspace(double xmin, double xmax, int N){
    double dx = (xmax - xmin) / (N-1);
    values = new double[N];
    for (int i = 0; i<N; i++)
        values[i] = xmin + dx * i;
}

linspace::~linspace()
{delete[] values;}

/****************************************
/ We need the constructor, linspace::linspace(double, double, int)
/ The constructor needs to create the values array (see, e.g., dep_vars::dep_vars(int)), then set the values of 
/ N, dx, and values
/
/ We need the desctructor, linspace::~linspace(). See, e.g., dep_vars::~dep_vars().
******************************************/