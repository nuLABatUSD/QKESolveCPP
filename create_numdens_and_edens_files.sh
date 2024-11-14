#!/usr/bin/bash

set -x

g++ number_density.cc QKE_methods.cc array_methods.cc matrices.cc thermodynamics.cc -o mon

./mon "T=1_QKESolveMPI_results.csv" "T=1_energydens.csv" "T=1_numberdens.csv"

./mon "not_equilibrium_QKESolveMPI_results.csv" "not_equilibrium_energydens.csv" "not_equilibrium_numberdens.csv"