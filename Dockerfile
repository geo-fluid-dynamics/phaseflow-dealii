FROM dealii/dealii:v8.5.0-gcc-mpi-fulldepscandi-debugrelease

RUN git clone https://github.com/alexanderzimmerman/phaseflow.git && \
    cd phaseflow && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make test && \
    cd ~
