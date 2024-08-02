#!/usr/bin/bash

set -x

mpic++ test2.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -std=c++11 -o wed

mpiexec -n 4 wed | tee 4coreoutput.txt

mpiexec -n 5 wed | tee 5coreoutput.txt

mpiexec -n 6 wed | tee 6coreoutput.txt

mpiexec -n 7 wed | tee 7coreoutput.txt

mpiexec -n 8 wed | tee 8coreoutput.txt

mpiexec -n 9 wed | tee 9coreoutput.txt

mpiexec -n 10 wed | tee 10coreoutput.txt

mpiexec -n 11 wed | tee 11coreoutput.txt

mpiexec -n 12 wed | tee 12coreoutput.txt