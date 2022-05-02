# Sliced-Optimal-Transport-Sampling

Source code of the sampler proposed in [Sliced Optimal Transport Sampling](https://perso.liris.cnrs.fr/lpaulin/Publications/paulin2020.html), Loïs Paulin, Nicolas Bonneel, David Coeurjolly, Jean-Claude Iehl, Antoine Webanck, Mathieu Desbrun, Victor Ostromoukhov, ACM Trans. on Graphics, SIGGRAPH 2020.

``` bibtex
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

This demo code generates unfiorm samples using the sliced optimal transport energy.


[![linux/macOS CI](https://github.com/loispaulin/Sliced-Optimal-Transport-Sampling/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/loispaulin/Sliced-Optimal-Transport-Sampling/actions/workflows/c-cpp.yml)

Dependencies:
=============
 + GSL
 + OpenMP (`brew install libomp` on macOS)

Code compilation:
=================

    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make sot


Code execution:
===================

    ./sot -n 1024 -d 2 -c -o test.dat

Generates 1 set of 1024 samples in dimension 2, stored in test.dat using default parameters

Detailed options:

```
%> sot -h
Help:
Option list:
	-o <OutputFileName> (optional): Specifies an output file in which points will be written.If unset standard output will be used
	-n <nbPoints> (default 1024): Specifies the number of points to generate
	-m <nbRealisations> (default 1): Specifies the number of generated pointsets
	-p <nbIteration> (default 4096): Specifies the number of batches in the optimization process
	--step <nbDirectionPerStep> (default 32): Specifies the number of slices per batch in the optimization process
	-s <seed> (default 133742): Specifies the random seed
	-d <dimension> (default 2): Specifies samples dimension
	-c (optional): If unset points will be given in the unit ball. Else they will be in the unit cube [0,1)^d
	--silent (optional): Cancels all outputs other than the points and errors

sot [-o <OutputFileName>] [-n <nbPoints>] [-m <nbRealisations>] [-p <nbIteration>][--step <nbDirectionPerStep>] [-s <seed>] [-d <dimension>] [-c]
``` 
