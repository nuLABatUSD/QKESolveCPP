#!/usr/bin/bash

set -x

mpic++ test4.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -std=c++11 -o wed

numprocs="64"
xmin="0."
xmax="10."
numlin="201"
numgl="5"
eta_e="0.01"
eta_mu="-0.01"
sin2theta="0.8"
deltamsquared="2.5e-15"
x_initial="0."
dx_initial="1.e12"
eta_relationship="same"
outputs="f_results_64_cores.csv"
times="times_64_cores.csv"

echo "#!/bin/bash" > cores_vs_eta_e_official_test_64_cores.sh
echo "# FILENAME: cores_vs_eta_e_test" >> cores_vs_eta_e_official_test_64_cores.sh 
echo "#SBATCH -A phy240216" >> cores_vs_eta_e_official_test_64_cores.sh 
echo "#SBATCH --nodes=1" >> cores_vs_eta_e_official_test_64_cores.sh 
echo "#SBATCH --ntasks=64" >> cores_vs_eta_e_official_test_64_cores.sh 
echo "#SBATCH -J qkesolvef" >> cores_vs_eta_e_official_test_64_cores.sh 
echo "#SBATCH -p shared" >> cores_vs_eta_e_official_test_64_cores.sh 
echo "#SBATCH --time=00:05:00" >> cores_vs_eta_e_official_test_64_cores.sh 
echo "#SBATCH --mail-user=lfowler@sandiego.edu" >> cores_vs_eta_e_official_test_64_cores.sh 
echo "#SBATCH--mail-type=all" >> cores_vs_eta_e_official_test_64_cores.sh

echo "module purge" >> cores_vs_eta_e_official_test_64_cores.sh 
echo "module load modtree/cpu" >> cores_vs_eta_e_official_test_64_cores.sh
echo "mpiexec -n $numprocs wed $xmin $xmax $numlin $numgl $eta_e $eta_mu $sin2theta $deltamsquared $x_initial $dx_initial $eta_relationship \"$outputs\" \"$times\"" >> cores_vs_eta_e_official_test_64_cores.sh