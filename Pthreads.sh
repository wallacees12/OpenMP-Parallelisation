#!/bin/bash

read -p "Enter number of threads: " THREADS
CC="gcc"
CFLAGS="-O3 -ftree-vectorize -march=native -ffast-math"
LDFLAGS="-lpthread"
NAME="galsim"
TARGET="./Pthreads/${NAME}"
SRCS="$TARGET.c" 

$CC $CFLAGS $SRCS -o $TARGET $LDFLAGS

time ./Pthreads/galsim 2000 ./input_data/ellipse_N_02000.gal 200 1e-5 0 $THREADS

