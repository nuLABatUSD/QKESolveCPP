#!/usr/bin/bash

set -x
rm wed

#input parameters will be xmin (1), xmax (2), numlin (3), numgl (4), eta_e (5), eta_mu (6), sin2theta (7), deltamsquared (8), N_step (9), dN (10), x_initial (11), x_final (12), dx_initial (13), verbose (14), output file name (15), density file name (16)

numprocs="4"
xmin="0."
xmax="10."
numlin="201"
numgl="5"
eta_e="0.01"
eta_mu="-0.01"
sin2theta="0.8"
deltamsquared="2.5e-15"
T="1"
N_step="1"
dN="1"
x_initial="0."
x_final="5.e15"
dx_initial="1.e12"
verbose="true"
output_file="T=1_QKESolveMPI_results.csv"
input_file="initialdensity.csv"

mpic++ static_main_QKESolveMPI.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -std=c++11 -o wed

mpiexec -n $numprocs wed $xmin $xmax $numlin $numgl $eta_e $eta_mu $sin2theta $deltamsquared $N_step $dN $x_initial $x_final $dx_initial $verbose $output_file $input_file