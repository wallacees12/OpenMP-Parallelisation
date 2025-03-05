#!/bin/bash

read -p "Enter number of threads: " THREADS
CC="gcc"
CFLAGS="-O3 -ftree-vectorize -march=native -ffast-math"
LDFLAGS="-lpthread -lm"
NAME="galsim_OpenMP"
DIR="Openmp/Joel"
TARGET="./$DIR/${NAME}"
SRCS="$TARGET.c" 
CMP="./compare_gal_files/compare_gal_files"
CMPSRC="$CMP.c"
CMPFILE="galsim_OpenMP.gal"

$CC $CFLAGS $SRCS -o $TARGET $LDFLAGS
$CC  $CMPSRC -o $CMP $LDFLAGS

$TARGET 2000 ./input_data/ellipse_N_02000.gal 200 1e-5 0 $THREADS

$CMP 2000 $CMPFILE ./ref_output_data/ellipse_N_02000_after200steps.gal

cp $CMPFILE $NAME.gal
mv $NAME.gal ./$DIR



