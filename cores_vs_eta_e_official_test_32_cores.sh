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


mpiexec -n 32 wed 0. 10. 201 5 0.01 -0.01 0.8 2.5e-15 0. 1.e12 same "results_32_cores.csv" "times_32_cores.csv"