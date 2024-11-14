#!/usr/bin/bash

g++ test3.cc alternative_integrals.cc matrices.cc thermodynamics.cc array_methods.cc QKE_methods.cc -o wed

./wed og_term1.csv og_term2.csv p3s.csv p4s.csv