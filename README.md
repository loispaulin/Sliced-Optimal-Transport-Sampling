# Sliced-Optimal-Transport-Sampling

Dependancies:
=============
 + GSL
 + OpenMP

Code compilation:
=================

    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make sot


Code toy execution:
===================

./sot -n 1024 -c -o test.dat

Generates 1 set of 1024 samples in test.dat using default parameters
