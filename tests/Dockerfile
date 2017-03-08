FROM dealii/dealii:v8.4.2-gcc-mpi-fulldepsmanual-release

MAINTAINER zimmerman.alex.g@gmail.com

RUN git clone https://github.com/alexanderzimmerman/nsb-pcm.git && \
    cd nsb-pcm && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    cd ~
