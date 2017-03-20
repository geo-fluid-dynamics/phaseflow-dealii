# nsb-pcm

The development of this code has just begun, starting from a copy of the author's convection-diffusion code at https://github.com/alexanderzimmerman/peclet

nsb-pcm (Navier-Stokes-Boussinesq with phase-change material) will solve the incompressible NSB equations coupled with an enthalpy-based energy equation allowing for melting and solidification of the phase-change material medium.

This is written in C++ and based on the deal.II finite element method library.

Author: Alexander G. Zimmerman <zimmerman@aices.rwth-aachen.de>

[![Build Status](https://travis-ci.org/alexanderzimmerman/nsb-pcm.svg?branch=master)](https://travis-ci.org/alexanderzimmerman/nsb-pcm) <- Continuous integration status (Click the button to go to Travis-CI)

# Run pre-built version on docker image

Pull the image from https://hub.docker.com/r/zimmerman/nsb-pcm/

```shell
docker pull zimmerman/nsb-pcm
```

# Build

This is currently being tested with the following builds of deal.II:
- deal.II v8.4.2 built by candi (https://github.com/koecher/candi) on Ubuntu 14.04
- deal.II v8.4.2 from docker image dealii/dealii:v8.4.2-gcc-mpi-fulldepsmanual-release (as shown in nsb-pcm/tests/Dockerfile)

To build:

```shell
git clone git@github.com:alexanderzimmerman/nsb-pcm.git

mkdir build

cd build

cmake ../nsb-pcm

make test
```
