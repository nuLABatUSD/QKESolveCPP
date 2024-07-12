#include "mpi.h"
#include <cmath>
#include <iostream>

using std::cout;
using std::cin;
using std::endl;

/* command line is 
mpic++ mpi_fun.cc -o wed
mpiexec -n 4 wed

*/

int main(int argc, char *argv[])
{

    int rows, cols;
    //need to initialize matrix type a, vector type b, and vector type c
    //also need vector buffer
    double** a = new double*[100];
    for(int k=0; k<100; k++){
        a[k] = new double[100];
    }
    double* b = new double[100];
    double* c = new double[100];
    double* buffer = new double[100];
    double ans;
    
    int myid, master, numprocs;
    MPI_Status* status;  
    int i, j, numsent, sender; 
    int anstype, row;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    master = 0;
    rows = 100;
    cols = 100;
    
    //note that master=0, this is the process performed only by the main program
    if (myid==master){
        //initialize a and b with random values
        for(i=0; i<cols; i++){
            b[i] = 1;
            for(j=0; j<rows; j++){
                a[i][j] = i;
            }
        }
        numsent = 0;
        //sending b, an object of length cols and type double, from master, to all other processors
        MPI_Bcast(b, cols, MPI_DOUBLE, master, MPI_COMM_WORLD);
        
        //send a row to each slave process; tag with row number (tag tells code where in c answer should be stored)
        for(i=0; i<std::min(numprocs-1,rows);i++){
            for(j=0; j<cols; j++){
                buffer[j] = a[i][j];
            }
            MPI_Send(buffer, cols, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
            numsent++;
        }
        
        for(i=0; i<rows;i++){
            MPI_Recv(&ans, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status);
            sender = status->MPI_SOURCE;
            anstype = status->MPI_TAG;
            c[anstype] = ans;
            if(numsent<=rows){
                for(j=0; j<cols; j++){
                    buffer[j] = a[numsent+1][j];
                }
                MPI_Send(buffer, cols, MPI_DOUBLE, sender, numsent+1, MPI_COMM_WORLD);
                numsent++;
                
            }
            else{
                MPI_Send(MPI_BOTTOM, 0, MPI_DOUBLE, sender, 0, MPI_COMM_WORLD);
            }
            
        }
        
        
    }
    
    //this is the process performed by all workers
    else{
        MPI_Bcast(b, cols, MPI_DOUBLE, master, MPI_COMM_WORLD);
        
        while(1){
            MPI_Recv(buffer, cols, MPI_DOUBLE, master, MPI_ANY_TAG, MPI_COMM_WORLD, status);
            if(status->MPI_TAG==0){
                break;
            }
            else{
                row = status->MPI_TAG;
                ans = 0;
                for(i=0; i<cols; i++){
                    ans += buffer[i]*b[i];
                }
                MPI_Send(&ans, 1, MPI_DOUBLE, master, row, MPI_COMM_WORLD);
            }
        }
        
    }
    
    MPI_Finalize();
    for(int k=0; k<100; k++){
        delete[] a[k];
    }
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] buffer;
    return 0;
    
    
}

/*

//finding pi program
    
    //n=user inputted value
    //myid=identifier for an individual processor
    //numprocs=number of total processors
    //mypi=individual processor contribution to pi estimate
    //pi=overall pi estimate
    //i,h,sum,x=all used in loops below
    
    
    
    
    int n, myid, numprocs, i;
    double PI25DT = 3.141592653589793238462643;
    double mypi, pi, h, sum, x;
    
    //beginning mpi program w number of processors as inputted
    MPI_Init(&argc, &argv);
    //sets number of processors
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    //gives an id number to each processor
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    //the while (1) works because it is always true, so always runs unless the break is hit
    while (1) {
        if (myid == 0){
            cout << "enter the number of intervals (0 quits): " << endl;
            cin >> n;
        }
        
        //sends value of n to all processors
        //MPI_Bcast(buffer, count, MPI datatype, root, MPI comm)
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (n==0)
            break;
        else{
            h = 1.0 / (double) n;
            sum = 0.0;
            //the i+= numprocs is what assigns jobs to processors (each processor gets every myid-th job
            for (i=myid+1; i<=n; i+=numprocs){
                x = h * ((double)i - 0.5);
                sum += (4.0 / (1.0 + x*x));
                
            }
            mypi = h * sum;
            
            MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if (myid ==0)
                cout << "pi is approximately. " << pi << ", Error is " << fabs(pi - PI25DT) << endl;
            
        }
        MPI_Finalize();
        return 0;
        
    }*/