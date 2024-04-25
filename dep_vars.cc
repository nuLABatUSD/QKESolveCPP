#include <iostream>

using std::cout;
using std::endl;
using std::ostream;

class dep_vars
{
    protected:
        int N;
    
    public:
        double* values;
    
        dep_vars(int size)
        {
            N = size;
            values = new double[N]();
        }
    
        dep_vars(double* copy_me, int size)
        {
            N = size;
            values = new double[N]();
            for (int i = 0; i < N; i++)
                values[i] = copy_me[i];
        }
    
        dep_vars(dep_vars* copy_me)
        {
            N = copy_me->length();
            values = new double[N]();
            for (int i = 0; i < N; i++)
                values[i] = copy_me->values[i];
        }
    
        ~dep_vars()
        {delete[] values;}
    
        int length()
        {return N;}
                 
    
        void print_all()
        {
            for (int i = 0; i < N; i++)
                cout << values[i] << endl;
        }
    
        void print(int N_top = 3, int N_bot = 1)
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
    
        
};

void f(double x, dep_vars* y, dep_vars* z)
{
    double* der = z->values;
    der[0] = 1.;
    der[1] = x;
    der[2] = 1/(1+x);
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
    
    delete input;
    delete output;
    return 0;
}