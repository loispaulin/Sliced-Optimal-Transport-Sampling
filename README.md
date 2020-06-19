# Sliced-Optimal-Transport-Sampling

Source code of the sampler proposed in "Sliced Optimal Transport Sampling", Loïs Paulin, Nicolas Bonneel, David Coeurjolly, Jean-Claude Iehl, Antoine Webanck, Mathieu Desbrun, Victor Ostromoukhov, ACM Trans. on Graphics, SIGGRAPH 2020.

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

```
@article{paulin2020,
      author = {Loïs Paulin and Nicolas Bonneel and David Coeurjolly and Jean-Claude Iehl and Antoine Webanck and Mathieu Desbrun and Victor Ostromoukhov},
      journal = {ACM Transactions on Graphics (Proceedings of SIGGRAPH)},
      month = {July},
      number = {4},
      title = {Sliced Optimal Transport Sampling},
      volume = {39},
      year = {2020}
}
```
