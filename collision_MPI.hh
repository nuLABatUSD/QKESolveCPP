#ifndef __COLLISION_MPI_HH__
#define __COLLISION_MPI_HH__

#include "arrays.hh"
#include "QKE_methods.hh"

class collision_MPI{
    protected:
        int myid, numprocs;
        
        linspace_and_gl* eps;
        nu_nu_collision** int_objects;
        
    public:
        collision_MPI(int, int, linspace_and_gl*);
        ~collision_MPI();
        
        void calculate_R(density*, dep_vars**);
        void collision_term(density*, dep_vars**, dep_vars**);
};

#endif