#!/usr/bin/bash

python3 create_initial_dens_from_file.py T=2_classical_not_equilibrium_QKESolveMPI_results.csv initialdensity.csv

g++ test3.cc alternative_integrals.cc matrices.cc thermodynamics.cc array_methods.cc QKE_methods.cc -o wed

./wed og_term1.csv og_term2.csv p3s.csv p4s.csv initialdensity.csv