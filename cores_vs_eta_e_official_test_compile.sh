#!/usr/bin/bash

set -x

mpic++ test4.cc QKESolveMPI.cc array_methods.cc QKE_methods.cc thermodynamics.cc matrices.cc QKESolve.cc -std=c++11 -o wed
