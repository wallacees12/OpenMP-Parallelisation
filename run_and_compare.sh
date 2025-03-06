#!/bin/bash

read -p "Enter number of threads: " THREADS
CC="gcc"
CFLAGS="-O3 -ftree-vectorize -march=native -ffast-math"
LDFLAGS="-lpthread -lm -fopenmp"
NAME="galsim_OpenMPv3"
DIR="Openmp/Joel"
TARGET="./$DIR/${NAME}"
SRCS="$TARGET.c" 
CMP="./compare_gal_files/compare_gal_files"
CMPSRC="$CMP.c"
CMPFILE="galsim_OpenMP.gal"

echo $TARGET

$CC $CFLAGS $SRCS -o $TARGET $LDFLAGS
$CC  $CMPSRC -o $CMP $LDFLAGS

$TARGET 3000 ./input_data/ellipse_N_03000.gal 100 1e-5 0 $THREADS

$CMP 3000 $CMPFILE ./ref_output_data/ellipse_N_03000_after100steps.gal

cp $CMPFILE $NAME.gal
mv $NAME.gal ./$DIR



