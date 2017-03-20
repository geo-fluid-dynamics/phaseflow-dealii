# nsb-pcm

The development of this code has just begun, starting from a copy of the author's convection-diffusion code at https://github.com/alexanderzimmerman/peclet

nsb-pcm (Navier-Stokes-Boussinesq with phase-change material) will solve the incompressible NSB equations coupled with an enthalpy-based energy equation allowing for melting and solidification of the phase-change material medium.

This is written in C++ and based on the deal.II finite element method library.

Author: Alexander G. Zimmerman <zimmerman@aices.rwth-aachen.de>

Doxygen generated HTML documentation: https://alexanderzimmerman.github.io/nsb-pcm/

[![Build Status](https://travis-ci.org/alexanderzimmerman/nsb-pcm.svg?branch=master)](https://travis-ci.org/alexanderzimmerman/nsb-pcm) <- Continuous integration status (Click the button to go to Travis-CI)

# For users:
## Run pre-built version on docker image

Pull the image from https://hub.docker.com/r/zimmerman/nsb-pcm/

    docker pull zimmerman/nsb-pcm

# For developers:
## Versions

This is currently being tested with the following builds of deal.II:
- deal.II v8.4.2 built by candi (https://github.com/koecher/candi) on Ubuntu 14.04
- deal.II v8.4.2 from docker image dealii/dealii:v8.4.2-gcc-mpi-fulldepsmanual-release (as shown in nsb-pcm/tests/Dockerfile)

## Build

    git clone git@github.com:alexanderzimmerman/nsb-pcm.git

    mkdir build

    cd build

    cmake ../nsb-pcm

    make test
    
## Documentation
The Doxygen generated HTML docs are hosted in the standard GitHub fashion on the gh-pages branch.

For initial set up with your local clone, follow the procedure found in https://github.com/m-a-d-n-e-s-s/madness/issues/104

Then whenever commiting to the master branch, update the gh-pages branch as follows:

    git push origin master

    doxygen

    cd doc/html

    git add *

    git commit -m "Refreshed HTML doc"

    git push origin gh-pages
