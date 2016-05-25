#!/bin/bash 
MAKEFLAGS="CC=mpicc" R CMD SHLIB -o mpihello_fromR.so mpihello_fromR.c --preclean
R -f mpihello.R
mpiexec -n 4 R -f mpihello.R