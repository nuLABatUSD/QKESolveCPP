#!/usr/bin/bash

set -x
rm wed

#input parameters will be xmin (1), xmax (2), numlin (3), numgl (4), eta_e (5), eta_mu (6), sin2theta (7), deltamsquared (8), N_step (9), dN (10), x_initial (11), x_final (12), dx_initial (13), verbose (14), output file name (15), density file name (16)

numprocs="128"
xmin="0."
xmax="20."
numlin="201"
numgl="5"
A="1"
B="8"
sin2theta="0"
deltamsquared="0"
T="16.0"
N_step="200"
dN="1"
x_initial="0."
x_final="5.e18"
dx_initial="1.e10"
verbose="true"
output_file="1e14_T=16_classical_not_equilibrium_QKESolveMPI_results.csv"
input_file="initialdensity.csv"
constants_output_file="1e14_T=16_classical_not_equilibrium_QKESolveMPI_constants.hh"

g++ initialize_dens.cc array_methods.cc matrices.cc thermodynamics.cc QKE_methods.cc -std=c++11 -o makedens
./makedens $xmin $xmax $numlin $numgl $T $A $B $input_file

echo "int numprocs = $numprocs;" > $constants_output_file
echo "double xmin = $xmin;" >> $constants_output_file
echo "double xmax = $xmax;" >> $constants_output_file
echo "int numlin = $numlin;" >> $constants_output_file
echo "int numgl = $numgl;" >> $constants_output_file
echo "int A = $A;" >> $constants_output_file
echo "int B = $B;" >> $constants_output_file
echo "double sin2theta = $sin2theta;" >> $constants_output_file
echo "double delamsquared = $deltamsquared;" >> $constants_output_file
echo "double T = $T;" >> $constants_output_file
echo "int N_step = $N_step;" >> $constants_output_file
echo "double dN = $dN;" >> $constants_output_file
echo "double x_initial = $x_initial;" >> $constants_output_file
echo "double x_final = $x_final;" >> $constants_output_file
echo "double dx_initial = $dx_initial;" >> $constants_output_file
echo "bool verbose = $verbose;" >> $constants_output_file
echo "std::string output_file_name = $output_file;" >> $constants_output_file

mpic++ static_main_QKESolveMPI.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -std=c++11 -o wed

echo "#!/bin/bash" > QKESolveMPI_execute.sh
echo "# FILENAME: 1MeV" >> QKESolveMPI_execute.sh 
echo "#SBATCH -A phy240216" >> QKESolveMPI_execute.sh 
echo "#SBATCH --nodes=1" >> QKESolveMPI_execute.sh 
echo "#SBATCH --ntasks=128" >> QKESolveMPI_execute.sh 
echo "#SBATCH -J qkesolvef" >> QKESolveMPI_execute.sh 
echo "#SBATCH -p wholenode" >> QKESolveMPI_execute.sh 
echo "#SBATCH --time=12:00:00" >> QKESolveMPI_execute.sh 
echo "#SBATCH --mail-user=lfowler@sandiego.edu" >> QKESolveMPI_execute.sh 
echo "#SBATCH--mail-type=all" >> QKESolveMPI_execute.sh

echo "module purge" >> QKESolveMPI_execute.sh 
echo "module load modtree/cpu" >> QKESolveMPI_execute.sh
echo "mpiexec -n $numprocs wed $xmin $xmax $numlin $numgl $sin2theta $deltamsquared $N_step $dN $x_initial $x_final $dx_initial $verbose \"$output_file\" \"$input_file\" " >> QKESolveMPI_execute.sh