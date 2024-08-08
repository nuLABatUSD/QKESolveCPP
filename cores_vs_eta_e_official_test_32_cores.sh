#!/bin/bash

# FILENAME: cores_vs_eta_e_test

#SBATCH -A phy240216
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH -J qkesolvef
#SBATCH -p shared
#SBATCH --time=00:05:00
#SBATCH --mail-user=lfowler@sandiego.edu
#SBATCH--mail-type=all

module purge
module load modtree/cpu


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
outputs="f_results_32_cores.csv"
times="times_32_cores.csv"

mpiexec -n 32 wed $xmin $xmax $numlin $numgl $eta_e $eta_mu $sin2theta $deltamsquared $x_initial $dx_initial $eta_relationship $ouputs $times
