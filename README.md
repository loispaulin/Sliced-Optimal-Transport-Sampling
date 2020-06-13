# Sliced-Optimal-Transport-Sampling

Dependancies:
=============
 + GSL
 + OpenMP (`brew install libomp` on macOS)

Code compilation:
=================

    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make sot


Code toy execution:
===================

    ./sot -n 1024 -d 2 -c -o test.dat

Generates 1 set of 1024 samples in dimension 2, stored in test.dat using default parameters
