#!/usr/bin/bash

set -x

mpic++ test4.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -std=c++11 -o wed

mpiexec -n 4 wed opposite "4coreoutputs.csv" "4coretimes.csv"

mpiexec -n 16 wed opposite "16coreoutputs.csv" "16coretimes.csv"

mpiexec -n 32 wed opposite "32coreoutputs.csv" "32coretimes.csv"

mpiexec -n 60 wed opposite "64coreoutputs.csv" "64coretimes.csv"