# phaseflow

The development of this code has just begun, starting from a copy of the author's convection-diffusion code at https://github.com/alexanderzimmerman/peclet

Phaseflow will solve the incompressible Navier-Stokes-Boussinesq equations coupled with an enthalpy-based energy equation allowing for melting and solidification of the phase-change material medium.

This is written in C++ and based on the deal.II finite element method library.

Author: Alexander G. Zimmerman <alex.g.zimmerman@gmail.com>

Doxygen generated HTML documentation: https://alexanderzimmerman.github.io/phaseflow/

[![Build Status](https://travis-ci.org/alexanderzimmerman/phaseflow.svg?branch=master)](https://travis-ci.org/alexanderzimmerman/phaseflow) (<b>Continuous integration status</b>; click the button to go to Travis-CI)

# For users:
## Run pre-built version on docker image
Get the free community edition of Docker here: https://www.docker.com/community-edition

Pull the image from https://hub.docker.com/r/zimmerman/phaseflow/ and run the container with docker

    docker run -ti zimmerman/phaseflow:latest
    
Or run the container with access to a shared folder (shared between the host and the container)

    docker run -ti -v $(pwd):/home/dealii/shared zimmerman/phaseflow:latest
    
If you plan to use this container repeatedly, then instead use this command to also give it a name

    docker run -ti -v $(pwd):/home/dealii/shared --name phaseflow zimmerman/phaseflow:latest

After exiting the container, you can start it again with

    docker start phaseflow
    
You can confirm that the container is running with

    docker ps
    
or list all containers (running or not) with

    docker ps -a

To enter a bash terminal inside of the running container

    docker start phaseflow
    
    docker exec -ti -u dealii phaseflow /bin/bash -l

# For developers:
## Versions

This is currently being tested with the following builds of deal.II:
- deal.II v8.5.0 from docker image dealii/dealii:v8.5.0-gcc-mpi-fulldepscandi-debugrelease

## Build

    git clone git@github.com:alexanderzimmerman/phaseflow.git

    mkdir build

    cd build

    cmake ../phaseflow

    make test
    
## Design notes
The Phaseflow class is implemented entirely with header files. This reduces the structural complexity of the code and can increase programming productivity, but it leads to longer compile times. A header-only approach would be impractical for the deal.II library itself; but in this small project's experience, the header-only approach is more than adequate. Most notably, this simplifies working with C++ templates.

## Documentation
The Doxygen generated HTML docs are hosted in the standard GitHub fashion on the gh-pages branch.

The procedure for keeping the HTML docs updated is rough. To initially create the gh-pages branch, we followed the outline in an [issue](https://github.com/m-a-d-n-e-s-s/madness/issues/104) from another repository, which uses the ideas from [here](http://rickfoosusa.blogspot.de/2011/10/howto-use-doxygen-with-github.html) and [here](https://gist.github.com/chrisjacob/825950). Since any of these links may break, in short the procedure was

    cd phaseflow
    
    mkdir doxygen_output
    
    mkdir doxygen_output/html
    
    cd doxygen_output/html
    
    git clone git@github.com:alexanderzimmerman/phaseflow.git .
    
    git checkout -b gh-pages
    
    git branch -d master
    
    git rm -r *

    git commit "Removed everything from gh-pages branch"
    
    cd ../..
    
    git checkout master
    
    doxygen
    
    cd doxygen_output/html
    
    git checkout gh-pages
    
    git add *
    
    git commit "Added all HTML documentation"
    
    git push origin gh-pages

You can skip many of those steps for initial set up with your local clone. Simply make the target documents directory and clone the gh-pages branch inside of it.

    cd phaseflow

    mkdir doxygen_output
    
    mkdir doxygen_output/html
    
    cd doxygen_output/html
    
    git clone git@github.com:alexanderzimmerman/phaseflow.git .
    
    git checkout gh-pages

Then whenever commiting to the master branch, update the gh-pages branch as follows:

    cd phaseflow

    doxygen

    cd doxygen_output/html

    git add *

    git commit -m "Refreshed HTML doc"

    git push origin gh-pages
