#include <iostream>
#include "constants.hh"

using std::cout;
using std::endl;

int main(){
    
    double*** M = new double**[4];
    for(int i =0; i<4; i++){
        M[i] = new double*[4];
    
        for(int j=0; j<4; j++){
            M[i][j] = new double[4];
        
            for(int k=0; k<4; k++){
               M[i][j][k] = i*j*k; 
            }
        }
    }
    
    cout << M[1][3][2];
    
    for(int i =0; i<4; i++){
        for(int j=0; j<4; j++){
           delete[] M[i][j]; 
        }
        delete[] M[i];
    }
    delete[] M;
    double today = 0;
    
    for(int i=6; i<today; i++){
        cout << "hooray";
        cout << i << endl;
    }
    
    return 0;
}
