#include <iostream>
#include "mpi.h"

using std::cout;
using std::endl;

double f(double, double);


int main(int argc, char *argv[]){
    int myid, numprocs;
    double x1;
    int y0_size, anstype, source;
    double* xvals = new double[6];
    double** yvals = new double*[6];
    double myans;
    MPI_Status* status;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    
    /* fifth order RK
    k1 = h*f(x,y)
    k2 = h*f(x+a2*h, y+b21*k1)
    k3 = h*f(x+a3*h, y+b31*k1+b32*k2)
    k4 = h*f(x+a4*h, y+b41*k1+b42*k2+b43*k3)
    k5 = h*f(x+a5*h, y+b51*k1+b52*k2+b53*k3+b54*k4)
    k6 = h*f(x+a6*h, y+b61*k1+b62*k2+b63*k3+b63*k4+b65*k5)
    y1 = y + c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5 + c6*k6
    
    */
    if(myid == 0){
        double h = 0.1;
        double a2 = 1/5;
        double a3 = 3/10;
        double a4 = 3/5;
        double a5 = 1;
        double a6 = 7/8;
        double b21 = 1/5;
        double b31 = 3/40;
        double b32 = 9/40;
        double b41 = 3/10;
        double b42 = -9/10;
        double b43 = 6/5;
        double b51 = -11/54;
        double b52 = 5/2;
        double b53 = -70/27;
        double b54 = 35/27;
        double b61 = 1631/55296;
        double b62 = 175/512;
        double b63 = 575/13824;
        double b64 = 44275/110592;
        double b65 = 253/4096;
        double c1 = 37/378;
        double c2 = 0;
        double c3 = 250/621;
        double c4 = 125/594;
        double c5 = 0;
        double c6 = 512/1771;

        y0_size = 3;
        x1 = 0.1;
        double* y0 = new double[y0_size];
        for(int i=0; i<y0_size; i++){
            y0[i] = i;
        }
        
        double* k1 = new double[y0_size];
        double* k2 = new double[y0_size];
        double* k3 = new double[y0_size];
        double* k4 = new double[y0_size];
        double* k5 = new double[y0_size];
        double* k6 = new double[y0_size];
        double* y1 = new double[y0_size];
        
        
        for(int i=0; i<6; i++){
            yvals[i] = new double[3];
        }
        
        xvals[0] = x1;
        xvals[1] = x1+a2*h;
        xvals[2] = x1+a3*h;
        xvals[3] = x1+a4*h;
        xvals[4] = x1+a5*h;
        xvals[5] = x1+a6*h;
        //yvals[0] = y0;
       
        for(int i=0; i<3; i++){
            yvals[0][i] = y0[i];
            yvals[1][i] = y0[i]+b21*k1[i];
            yvals[2][i] = y0[i]+b31*k1[i]+b32*k2[i];
            yvals[3][i] = y0[i]+b41*k1[i]+b42*k2[i]+b43*k3[i];
            yvals[4][i] = y0[i]+b51*k1[i]+b52*k2[i]+b53*k3[i]+b54*k4[i];
            yvals[5][i] = y0[i]+b61*k1[i]+b62*k2[i]+b63*k3[i]+b64*k4[i]+b65*k5[i];
        }

        double** kvals = new double*[6];
        kvals[0] = k1;
        kvals[1] = k2;
        kvals[2] = k3;
        kvals[3] = k4;
        kvals[4] = k5;
        kvals[5] = k6;
        
        
        
        for(int i=0; i<6; i++){
            for(int j=1; j<y0_size; j++){
                MPI_Recv(&myans, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status);
                source = status->MPI_SOURCE;
                anstype = status->MPI_TAG;
                kvals[i][anstype] = h*myans;
            }
        }
        
        for(int i=0; i<y0_size; i++){
            y1[i] = yvals[0][i] + c1*k1[i] + c2*k2[i] + c3*k3[i] + c4*k4[i] + c5*k5[i] + c6*k6[i]; 
        }
        
        
        cout << "This is y1: " << endl;
        for(int i=0; i<y0_size; i++){
            cout << y1[i] << endl;
        }
        
  
        delete[] y0;
        delete[] y1;
        delete[] xvals;
        
        for(int i=0; i<6; i++){
            delete[] kvals[i];
        }
        delete[] kvals;
    
     
    }
    
    
    
    for(int i=0; i<6; i++){
        //everybody gets the y for that step
        MPI_Bcast(&xvals[i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(yvals[i], y0_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if(myid != 0){
            for(int j=0; j<y0_size; j+=numprocs){
                myans = f(xvals[i], yvals[i][j]);
                MPI_Send(&myans, 1, MPI_DOUBLE, 0, j, MPI_COMM_WORLD);
                
            }
        }
        /*
        if(myid == 0){
            for(int j=1; j<y0_size; j++){
                MPI_Recv(&myans, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status);
                anstype = status->MPI_TAG;
                kvals[i][anstype] = h*myans;
            }
        }*/
    }
    
    /*if(myid == 0){
        for(int i=0; i<y0_size; i++){
            y1[i] = yvals[0][i] + c1*k1[i] + c2*k2[i] + c3*k3[i] + c4*k4[i] + c5*k5[i] + c6*k6[i]; 
        }
        
        
        cout << "This is y1: " << endl;
        for(int i=0; i<y0_size; i++){
            cout << y1[i] << endl;
        }
        
  
        delete[] y0;
        delete[] y1;
        delete[] xvals;
        
        for(int i=0; i<6; i++){
            delete[] kvals[i];
        }
        delete[] kvals;
    }*/
    
    
    
    
    for(int i=0; i<6; i++){
        delete[] yvals[i];
    }
    delete[] yvals;
    delete[] xvals;
    MPI_Finalize();
    
}


double f(double x, double y){
    return x*y*0.3;
    
}