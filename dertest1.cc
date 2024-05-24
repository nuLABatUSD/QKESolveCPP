// import statements
#include <iostream>
#include <cmath>
#include <algorithm>
#include "constants.hh" // so now I can just write the names of variables in constants.hh
#include "CashKarp_vals.hh"
#include "arrays.hh"
#include <string>
#include <fstream>

using std::cout;
using std::endl;
using std::ofstream;

void f(double, dep_vars*, dep_vars*);
void RKCash_Karp(double, dep_vars*, double, double*, dep_vars*, dep_vars*);
bool step_accept(dep_vars*, dep_vars*, dep_vars*, double, double*);
bool RKCK_step(double, dep_vars*, double, double*, dep_vars*, double*);
bool ODEOneRun(double, dep_vars*, double, int, int, double, double*, dep_vars*, double*, const std::string& file_name);





int main(int argc, char** argv) 
{     
    int N;
    if (argc == 5) 
    {
        N = 4;
        //cout << "Using der2 " << endl;
        
    } else { // NOTE: for some reason you have to put the } else { on the same line
        //cout << "Using der1 " << endl;
        N = 3;
    }
    
    //initialize variables
    double x = 0;
    double x_stepped = 0;
    double dx = 0.01;
    
    dep_vars* y = new dep_vars(N);
    dep_vars* y_5th = new dep_vars(N);
    dep_vars* y_4th = new dep_vars(N);
    dep_vars* der = new dep_vars(N);
    dep_vars* y_next = new dep_vars(N);
    double* x_next = new double[N]();
    double* dx_next = new double[N]();
    dep_vars* y5 = new dep_vars(N);
    dep_vars* y4 = new dep_vars(N);
    double* dx_new = new double[N]();
    
    double x0 = 0.0;
    double dx0 = 0.01;
    int N_step = 1000;
    int dN = 10;
    double x_final = 10.0;
    double x_array; 
    dep_vars y_array(N);
    double dx_array; 
    std::string file_name = "output.csv";
    
    if (argc == 5)
    {       
       for (int i = 0; i<N; i++)
       {
           y -> set_value(i, atof(argv[i+1]));
       }
     }
    
    f(x, y, der);
    
    RKCash_Karp(x, y, dx, &x_stepped, y_5th, y_4th); // the & delivers the memory address , a pointer    

    RKCK_step(x, y, dx, x_next, y_next, dx_next);
    
    
    
    
// ODEOneRun make data
    ODEOneRun(x0, y, dx0, N_step, dN, x_final, &x_array, &y_array, &dx_array, file_name);
    //inputting a double means cant change this value in the function, vs messing around with y actually changes things when we get back to main (because its a pointer),,,, outputs are always pointers, becuase we want to change memory so function can do somethi
    // & "reference to" instead of giving number that x_array is equal to its giving us a reference or pointer to x_array... the & takes the thing and turns it into a pointer .. sends the memory address to the function so that the function can go at it with (when i have a data thing and want to create a pointer to the data thing, used when we call functions
    
    // delete them
    delete y;
    delete der;
    delete y_5th;
    delete y_4th;
    delete y_next;
    delete[] x_next;
    delete[] dx_next;
    delete y5;
    delete y4;
    delete[] dx_new;

    return 0; 
}














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











/***************************
double* x_stepped is a pointer to the memory spot of a single double precision number
double* y_5th and 4th are pointers to arrays of doubles each with N elements
***************************/
void RKCash_Karp(double x, dep_vars* y, double dx, double* x_stepped, dep_vars* y_5th, dep_vars* y_4th)
{
    //int N;
    int N = y->length();
    //to use k1 - need to delcaire it as an array of N doubles and allocate memory (remember to delete after)
    dep_vars* k1 = new dep_vars(N);
    dep_vars* k2 = new dep_vars(N);
    dep_vars* k3 = new dep_vars(N);
    dep_vars* k4 = new dep_vars(N);
    dep_vars* k5 = new dep_vars(N);
    dep_vars* k6 = new dep_vars(N);
    
    dep_vars* z2 = new dep_vars(N);
    dep_vars* z3 = new dep_vars(N);
    dep_vars* z4 = new dep_vars(N);
    dep_vars* z5 = new dep_vars(N);
    dep_vars* z6 = new dep_vars(N); //inputs to get k2-k6
    
    // k1 = dx * f(x, y)
    f(x, y, k1);
    k1 -> multiply_by(dx);  //k1 = dx * f(x,y)
  
    // k2 = dx * f(x + a2*dx, y + b21*k1)
    z2 -> copy(y);           //z2 = y
    z2 -> add_to(b21, k1);      //z2 = y + b21*k1
    f(x + a2*dx, z2, k2);          //k2 = f(x+a2*dx, z2)
    k2 -> multiply_by(dx);     //dx*f(..)
    
    // k3 = dx * f(x + a3*dx, y + b31*k1 + b32*k2)
    z3 -> copy(y);           //z3 = y
    z3 -> add_to(b31, k1); //z3 = y + b31*k1
    z3 -> add_to(b32, k2);
    f(x + a3*dx, z3, k3);         // k3 = f(x + a3*dx, z3)
    k3 -> multiply_by(dx);  // k3 = dx*f(x + a3*dx, z3)
 
    // k4 = dx * f(x + a4*dx, y + b41*k1 + b42*k2 +b43*k3)
    z4 -> copy(y);           //z4 = y
    z4 -> add_to(b41, k1);  //z4 = y + b41*k1
    z4 -> add_to(b42, k2); //z4 = y + b41*k1 + b42*k2
    z4 -> add_to(b43, k3); //z4 = y + b41*k1 + b42*k2 + b43*k3
    f(x + a4*dx, z4, k4);         // k4 = f(x + a4*dx, z4)
    k4 -> multiply_by(dx);
        
    // k5 = dx * f(x + a5*dx, y + b51*k1 + b52*k2 + b53*k3 + b54*k4)
    z5 -> copy(y);           //z5 = y
    z5 -> add_to(b51, k1);      //z5 = y + b51*k1
    z5 -> add_to(b52, k2);      //z5 = y + b51*k1 + b52*k2
    z5 -> add_to(b53, k3);      //z5 = y + b51*k1 + b52*k2 + b53*k3
    z5 -> add_to(b54, k4);      //z5 = y + b51*k1 + b52*k2 + b53*k3 + b54*k4    
    f(x + a5*dx, z5, k5);         // k5 = f(x + a5*dx, z5)
    k5 -> multiply_by(dx);
    
    // k6 = dx * f(x + a6*dx, y + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5)
    z6 -> copy(y);           //z6 = y
    z6 -> add_to(b61, k1);      //z6 = y + b61*k1
    z6 -> add_to(b62, k2);      //z6 = y + b61*k1 + b62*k2
    z6 -> add_to(b63, k3);      //z6 = y + b61*k1 + b62*k2 + b63*k3
    z6 -> add_to(b64, k4);      //z6 = y + b61*k1 + b62*k2 + b63*k3 + b64*k4
    z6 -> add_to(b65, k5);      //z6 = y + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5 
    f(x + a6*dx, z6, k6);         // k6 = f(x + a6*dx, z6)
    k6 -> multiply_by(dx);
     
    //y_5th = y + c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5 + c6*k6
    y_5th -> copy(y); //y_5th = y
    y_5th -> add_to(c1, k1);
    y_5th -> add_to(c2, k2);
    y_5th -> add_to(c3, k3);
    y_5th -> add_to(c4, k4);
    y_5th -> add_to(c5, k5);
    y_5th -> add_to(c6, k6);
            
    // y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3 + cstar4*k4 + cstar5*k5 + cstar6*k6
    y_4th -> copy(y); //y_4th = y           
    y_4th -> add_to(cstar1, k1); //y_4th = y + cstar1*k1
    y_4th -> add_to(cstar2, k2); //y_4th = y + cstar1*k1 + cstar2*k2
    y_4th -> add_to(cstar3, k3); //y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3
    y_4th -> add_to(cstar4, k4); //y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3 + cstar4*k4
    y_4th -> add_to(cstar5, k5); //y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3 + cstar4*k4 + cstar5*k5
    y_4th -> add_to(cstar6, k6); //y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3 + cstar4*k4 + cstar5*k5 + cstar6*k6
    
    // x_stepped = x + dx
    *x_stepped = x + dx;
    
         
    delete k1;
    delete k2;
    delete k3;
    delete k4;
    delete k5;
    delete k6;
   
    delete z2;
    delete z3;
    delete z4;
    delete z5;
    delete z6;

    return;
}


















double eps = 1e-8; 
double TINY = 1e-40;
double Safety = 0.9;



bool step_accept(dep_vars* y, dep_vars* y5, dep_vars* y4, double dx, double* dx_new)
{
    int N = y->length();
    
    double dsm = 0;
    double delta1 = 0;
    double delta0 = 0;

    for (int i = 0; i<N; i++)
    { 
        delta1 = abs(y5 -> get_value(i) - y4 -> get_value(i));
        delta0 = eps*(abs(y -> get_value(i)) + abs(y5 -> get_value(i) - y -> get_value(i))) + TINY;
        
        if (delta1/delta0 > dsm)
        { 
            dsm = delta1/delta0;
            
         }
     }
      
    if (dsm == 0)
    {
        *dx_new = 5 * dx;
        //cout<< "TRUE (dsm == 0) dx_new = " << *dx_new << endl;
        return true;
    } 
    else if (dsm < 1){
        *dx_new = Safety * dx * pow(dsm, -0.2);
        *dx_new = std::min(5.0 * dx, *dx_new); 
        //cout<< "TRUE (dsm < 1) dx_new = " << *dx_new << endl;
        return true;
    }
    else{
        *dx_new = Safety * dx * pow(dsm, -0.25);
        //cout<< "FALSE dx_new = " << *dx_new << endl;
        return false;
    }
    
    
}















bool RKCK_step(double x, dep_vars* y, double dx, double* x_next, dep_vars* y_next, double* dx_next)
{
    double dx_try = dx;
    int N = y->length();
    dep_vars* y5 = new dep_vars(N); //???
    dep_vars* y4 = new dep_vars(N);
    
    for (int i = 0; i<10; i++)
        
    { 
        RKCash_Karp(x, y, dx_try, x_next, y5, y4);
        if (step_accept(y, y5, y4, dx_try, dx_next))
        {
            y_next -> copy(y5);
            return true;
            break;
        } 
        else {
           dx_try = *dx_next; 
        }
        
    }
    
    cout << "ERROR:  10 iterations without acceptable step" << endl;
    
    delete y5;
    delete y4;
    
    return false;
    
    
}


/*********************************************************************************
***********************************************************************************

************************************************
HOW TO: WRITE A STRING AS AN INPUT TO A FUNCTION:
************************************************

#include <string>

// Function that takes a string as an argument
void printMessage(const std::string& message) {
    std::cout << "The message is: " << message << std::endl;
}

int main() {
    // Define a string
    std::string myMessage = "Hello, World!";

    // Call the function with the string as an argument
    printMessage(myMessage);

    return 0;
}

************************************************
********* EXAMPLE FROM MAIN.CC
************************************************
int main(){
    
    ofstream file("energy_density_and_pressure_for_T.csv");
    
    file << "T, density, pressure" << endl;
    
    for (int i=1; i<101; i++){
        double dens = 0.;
        double press =0.;
        energy_and_pressure(_electron_mass_, i, &dens, &press);
        
        file << i << ", " << dens << ", " << press << endl;
    }
    file.close();
    
}


************************************************
HOW TO: OPEN A FILE WITH A STRING NAME:
************************************************


#include <iostream>
#include <fstream>
#include <string>

int main() {
    std::string filename = "example.txt";
    
    // Create and open a text file for writing
    std::ofstream outFile(filename);
    if (outFile.is_open()) {
        outFile << "Hello, this is a test file.\n";
        outFile << "Writing some data to the file.\n";
        outFile.close(); // Close the file
    } else {
        std::cerr << "Unable to open file for writing: " << filename << std::endl;
        return 1;
    }
    
    // Open the same file for reading
    std::ifstream inFile(filename);
    if (inFile.is_open()) {
        std::string line;
        while (getline(inFile, line)) {
            std::cout << line << std::endl; // Print each line to the console
        }
        inFile.close(); // Close the file
    } else {
        std::cerr << "Unable to open file for reading: " << filename << std::endl;
        return 1;
    }
    
    return 0;
}
**********************************************************************************
********************************************************************************* 
*/
          

bool ODEOneRun(double x0, dep_vars* y0, double dx0, int N_step, int dN, double x_final, double* x, dep_vars* y, double* dx, const std::string& file_name) 
{
    // Set x, y, dx to initial values
    int N = y -> length();
    *x = x0;
    y -> copy(y0);
    *dx = dx0;

    // Declare for RKCK_step
    double* x_next = new double; 
    dep_vars* y_next = new dep_vars(N);
    double* dx_next = new double; 
    
    ofstream file(file_name);
    file << "x, dx, y1, y2, y3, y4" << endl;
    
    //print to the file - initial values
    file << *x << ", " << *dx << ", ";
    for (int k = 0; k < N; k++) 
    {
        //file << y -> get_value(k) << " ";
        file << y->get_value(k);
         if (k < N - 1) file << ", "; // had to make sure there was no comma at the end of the y array for the np.arary to work when loading the output data
    }
    file << endl;
    
    
    for (int i = 0; i < N_step; i++) 
    {
        for (int j = 0; j < dN; j++) 
        {
           // cout << "ith step: " << i << ", jth step: " << j << endl;

            if (*x + *dx > x_final) 
            {
                *dx = x_final - *x;
            }
            
            if (RKCK_step(*x, y, *dx, x_next, y_next, dx_next)) 
            {
                // Update x, y, dx with the results from the RKCK step
               // cout << "Before update... x =  " << *x << "and dx = " << *dx << endl;
              
                *x = *x_next;
                y->copy(y_next);
                *dx = *dx_next;
               
               // cout << "After update... x = " << *x << "and dx = " << *dx << endl;

               // Write the updated values to the file here if only want ALL vals
                
            } 
            else 
            {
                delete x_next;
                delete y_next;
                delete dx_next;
                file.close();
                return false;
            }

            if (*x == x_final) 
            {
                cout << "Reached x_final" << endl;
                // Write the updated values to the file here if only want final vals
                file << *x << ", " << *dx << ", ";
                
                for (int k = 0; k < N; k++) 
                {
                    //file << y -> get_value(k) << " ";
                     file << y->get_value(k);
                     if (k < N - 1) file << ", ";
                }
                file << endl;
                delete x_next;
                delete y_next;
                delete dx_next;
                file.close();
                return true;
            }
        }
        
        
        file << *x_next << ", " << *dx_next << ", ";
        for (int k = 0; k < N; k++) 
        {
             //file << y -> get_value(k) << " ";
             file << y->get_value(k);
             if (k < N - 1) file << ", ";
         }
        file << std::endl;
    }
    
    delete x_next;
    delete y_next;
    delete dx_next; 
    file.close();
    return true;
}

