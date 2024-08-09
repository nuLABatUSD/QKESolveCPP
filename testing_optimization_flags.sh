#!/usr/bin/bash

set -x

mpic++ test3.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -std=c++11 -o wed

mpiexec -n 8 wed | tee outputo0.txt

mpic++ test3.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -o1 -std=c++11 -o wed

mpiexec -n 8 wed | tee outputo1.txt

mpic++ test3.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -o2 -std=c++11 -o wed

mpiexec -n 8 wed | tee outputo2.txt

mpic++ test3.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -o3 -std=c++11 -o wed

mpiexec -n 8 wed | tee outputo3.txt