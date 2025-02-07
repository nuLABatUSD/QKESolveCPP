#!/usr/bin/bash

set -x

g++ number_density.cc QKE_methods.cc array_methods.cc matrices.cc thermodynamics.cc alternative_integrals.cc -std=c++2a -o mon

./mon "initialdensity.csv" "energydens.csv" "numberdens.csv"