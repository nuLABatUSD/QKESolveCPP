#!/usr/bin/bash

set -x

mpic++ test2.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -std=c++11 -o wed

mpiexec -n 4 wed | tee -a output.txt

mpiexec -n 5 wed | tee -a output.txt

mpiexec -n 6 wed | tee -a output.txt

mpiexec -n 7 wed | tee -a output.txt

mpiexec -n 8 wed | tee -a output.txt

mpiexec -n 9 wed | tee -a output.txt

mpiexec -n 10 wed | tee -a output.txt

mpiexec -n 11 wed | tee -a output.txt

mpiexec -n 12 wed | tee -a output.txt