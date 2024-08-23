#!/usr/bin/bash

set -x
rm wed

numprocs="8"
xmin="0."
xmax="10."
numlin="201"
numgl="5"
sin2theta="0"
deltamsquared="0"
eta_e="0.01"
eta_mu="-0.01"
T="2.0"
N_step="2"
dN="1"
x_initial="0."
x_final="5000000000000000000"
dx_initial="158983679795032800"
verbose="true"
previous_file="1e15_T=2_classical_not_equilibrium_QKESolveMPI_results.csv"
output_file="1e17_T=2_classical_not_equilibrium_QKESolveMPI_results.csv"
input_file="initialdensity.csv"

python3 create_initial_dens_from_file.py $previous_file $input_file

mpic++ static_main_QKESolveMPI.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -std=c++11 -o wed

mpiexec -n $numprocs wed $xmin $xmax $numlin $numgl $sin2theta $deltamsquared $N_step $dN $x_initial $x_final $dx_initial $verbose $output_file $input_file