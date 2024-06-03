#include <complex>

using std::complex;

class complex_three_vector{
    private:
        complex<double>* values;
        
    public:
    complex_three_vector(int Nv=3);
    complex_three_vector(complex<double>, complex<double>, complex<double>);
    complex_three_vector(complex<double>);
    complex_three_vector(complex_three_vector*);
    
    void print_all();
    complex<double> get_value(int);

    void add(complex_three_vector*, complex_three_vector*);
    complex<double> dot_with(complex_three_vector*);
    complex<double> magnitude_squared();
    complex<double> magnitude();
    void set_cross_product(complex_three_vector*, complex_three_vector*);
    
    ~complex_three_vector();
    
};
/*
class matrix{
    protected:
        complex A0;
        complex_three_vector A;
    
    public:
    
    
};
*/