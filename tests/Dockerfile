FROM dealii/dealii:v8.5.0-gcc-mpi-fulldepscandi-debugrelease

RUN git clone https://github.com/alexanderzimmerman/nsb-pcm.git && \
    cd nsb-pcm && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    cd ~
