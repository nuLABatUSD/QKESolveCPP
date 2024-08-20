#!/usr/bin/bash

set -x
rm wed

#input parameters will be xmin (1), xmax (2), numlin (3), numgl (4), eta_e (5), eta_mu (6), sin2theta (7), deltamsquared (8), N_step (9), dN (10), x_initial (11), x_final (12), dx_initial (13), verbose (14), output file name (15), density file name (16)

numprocs="128"
xmin="0."
xmax="10."
numlin="201"
numgl="5"
A="1"
B="8"
sin2theta="0"
deltamsquared="0"
T="2.0"
N_step="200"
dN="1"
x_initial="0."
x_final="5.e15"
dx_initial="1.e12"
verbose="true"
output_file="T=2_classical_not_equilibrium_QKESolveMPI_results.csv"
input_file="initialdensity.csv"

g++ initialize_dens.cc array_methods.cc matrices.cc thermodynamics.cc QKE_methods.cc -std=c++11 -o makedens
./makedens $xmin $xmax $numlin $numgl $T $A $B $input_file


mpic++ static_main_QKESolveMPI.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -std=c++11 -o wed

echo "#!/bin/bash" > QKESolveMPI_execute.sh
echo "# FILENAME: 1MeV" >> QKESolveMPI_execute.sh 
echo "#SBATCH -A phy240216" >> QKESolveMPI_execute.sh 
echo "#SBATCH --nodes=1" >> QKESolveMPI_execute.sh 
echo "#SBATCH --ntasks=128" >> QKESolveMPI_execute.sh 
echo "#SBATCH -J qkesolvef" >> QKESolveMPI_execute.sh 
echo "#SBATCH -p wholenode" >> QKESolveMPI_execute.sh 
echo "#SBATCH --time=02:00:00" >> QKESolveMPI_execute.sh 
echo "#SBATCH --mail-user=lfowler@sandiego.edu" >> QKESolveMPI_execute.sh 
echo "#SBATCH--mail-type=all" >> QKESolveMPI_execute.sh

echo "module purge" >> QKESolveMPI_execute.sh 
echo "module load modtree/cpu" >> QKESolveMPI_execute.sh
echo "mpiexec -n $numprocs wed $xmin $xmax $numlin $numgl $sin2theta $deltamsquared $N_step $dN $x_initial $x_final $dx_initial $verbose \"$output_file\" \"$input_file\" " >> QKESolveMPI_execute.sh