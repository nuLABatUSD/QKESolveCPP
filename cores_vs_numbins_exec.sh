#!/usr/bin/bash

set -x

mpic++ test2.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -std=c++11 -o wed

mpiexec -n 4 wed | tee -a output.txt

mpiexec -n 3 wed | tee -a output.txt