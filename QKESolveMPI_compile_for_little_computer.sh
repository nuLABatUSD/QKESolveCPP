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
T="0.25"
N_step="2"
dN="1"
x_initial="0."
x_final="5.e15"
dx_initial="1.e12"
verbose="true"
output_file="electron_and_neutrino_together_QKESolveMPI_results.csv"
input_file="initialdensity.csv"

g++ initialize_thermal_dens.cc array_methods.cc matrices.cc thermodynamics.cc QKE_methods.cc -std=c++11 -o makedens
./makedens $xmin $xmax $numlin $numgl $T $eta_e $eta_mu $input_file

mpic++ static_main_QKESolveMPI.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -std=c++11 -o wed

mpiexec -n $numprocs wed $xmin $xmax $numlin $numgl $sin2theta $deltamsquared $N_step $dN $x_initial $x_final $dx_initial $verbose $output_file $input_file