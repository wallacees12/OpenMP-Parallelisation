#!/bin/bash

read -p "Enter number of threads: " THREADS
CC="gcc"
CFLAGS="-O3 -ftree-vectorize -march=native -ffast-math"
LDFLAGS="-lpthread -lm"
NAME="galsim"
TARGET="./Pthreads/${NAME}"
SRCS="$TARGET.c" 
CMP="./compare_gal_files/compare_gal_files"
CMPSRC="$CMP.c"

$CC $CFLAGS $SRCS -o $TARGET $LDFLAGS
$CC  $CMPSRC -o $CMP $LDFLAGS

time $TARGET 2000 ./input_data/ellipse_N_02000.gal 200 1e-5 0 $THREADS

#./$CMP result.gal ./

