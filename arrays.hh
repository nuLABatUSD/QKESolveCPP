#ifndef _ARRAYS_HH_
#define _ARRAYS_HH_

class dep_vars
{
    protected:
        int N;
        double* values;
    
    public:
        /******************************
        /  Constructors:
        /  Create an object with:
        /  (int): specify number of dependent variables, initialize them to zero
        /  (double*, int): include a pointer to an array and the length, initialize
        /  (dep_vars*): copy another dep_vars object
        *******************************/
        dep_vars(int);
        dep_vars(double*, int);
        dep_vars(dep_vars*);
        ~dep_vars();

        /*****************************
        /  Return protected values: length and individual terms values[i]
        /  Protect the array of values by returning numbers, not a pointer to the array
        *****************************/
        int length();
        double get_value(int);

        /*****************************
        /  These methods actually change the values array. Use with caution.
        /  set_value(int i, double v) sets values[i] = v
        /  copy(dep_vars*) copies the elements of another dep_vars object
        /  multiply_by(double c) takes values ==> c * values
        /  add_to(double c, dep_vars* z) takes values ==> values + c * z
        ******************************/

        void set_value(int, double);
        void copy(dep_vars*);
        void multiply_by(double);
        void add_to(double, dep_vars*);

        /*****************************
        /  Methods to print to stdout
        /  print_all() prints every term on a new line (can be very long...)
        /  print(int N_top, int N_bottom) prints the first N_top terms (each on their own line) followed by "..." and then N_bottom terms.
        /  default is to print top 3, then final term... e.g., print() creates this behavior
        ******************************/

        void print_all();
        void print(int, int);
            
};



class three_vector : public dep_vars
{
    public:
    /*****************************
    / Constructors:
    / three_vector() - initializes to zeros
    / three_vector(double, double, double) - uses three values to initialize x, y, z
    / three_vector(double*) - uses the first three values of the array to initialize x, y, z
    /
    / Does not need its own destructor. Just uses dep_vars destructor.
    *****************************/
    three_vector();
    three_vector(double, double, double);
    three_vector(double*);
    three_vector(three_vector*);

    /******************************
    / dot_with(three_vector*), magnitude_squared(), and magnitude() return a double. They represent dot product 
    / of the vector with a second vector, itself, and the square root thereof, respectively.
    / set_cross_product(three_vector*, three_vector*) actually overwrites the components of this three_vector. Use with care.
    ******************************/

    double dot_with(three_vector*);
    double magnitude_squared();
    double magnitude();
    void set_cross_product(three_vector*, three_vector*);

};

class density : public dep_vars
{
    protected:
    int N_bins;
    
    public:
    density(int);
    density(linspace*, double, double);
    
    int num_bins();
    void p_vector(int, bool, three_vector*);
    void p0_p(int, bool, three_vector*);

};

struct dummy_vars{
    int N;
    double* values;
    double* dx;
    
    dummy_vars(int);
    void print_all();
    ~dummy_vars();
};
    
struct linspace : public dummy_vars
{
    linspace(double, double, int);
};

#endif