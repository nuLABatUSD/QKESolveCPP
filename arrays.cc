#include "arrays.hh"

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

/****************************************
/ We need functions
/ double three_vector::magnitude()
/ void three_vector::set_cross_product(three_vector* A, three_vector* B)
/
/ The latter sets the components of the three-vector (e.g., values[0], values[1], values[2]) with the components of AxB
/ Note that this over-writes whatever info lives within the three_vector object already.
****************************************/

/****************************************
/ We need the constructor, linspace::linspace(double, double, int)
/ The constructor needs to create the values array (see, e.g., dep_vars::dep_vars(int)), then set the values of 
/ N, dx, and values
/
/ We need the desctructor, linspace::~linspace(). See, e.g., dep_vars::~dep_vars().
******************************************/