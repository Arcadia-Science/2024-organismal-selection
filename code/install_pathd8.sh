#!/bin/bash
# Install PATHd8
mkdir tmpdir && cd tmpdir
wget https://www2.math.su.se/PATHd8/PATHd8.zip
unzip PATHd8.zip
cc PATHd8.c -O3 -lm -o PATHd8
mv PATHd8 $CONDA_PREFIX/bin
rm -r ./*
