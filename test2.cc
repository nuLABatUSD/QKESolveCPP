#include <iostream>
#include <time.h>

using namespace std; 
using std::cout;
using std::endl;

bool is_prime(int);

int main(int argc, char *argv[])
{
    clock_t t;
    t = clock();
    
    int n;
    cout << "enter the max number: " << endl;
    cin >> n;
   
    int mycount = 0;
    for (int i=2; i<=n; i++){
        //check all numbers it is responsible for then send those back to main
        if(is_prime(i)){
            mycount++;
        }
    }
    
    cout << "There are " << mycount << " primes less than " << n << "." << endl;
   
    t = clock() - t;
    cout << "Time it takes not parallelized " << t << " clock ticks" << endl;
    return 0;
    
}


bool is_prime(int p){
    bool is_prime = true;
    for(int i=2; i<p; i++){
        if(p%i == 0){
            is_prime = false;
            break;
        }
    }
    
    return is_prime;
    
}