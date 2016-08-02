FROM dealii/dealii:manual-debug

MAINTAINER luca.heltai@gmail.com

USER root

RUN apt-get update && apt-get -yq install \
    assimp-utils \
    libassimp-dev

USER $USER

#sundials
SUNDIALS_VERSION 2.6.2
RUN wget https://github.com/luca-heltai/sundials/archive/v$SUNDIALS_VERSION.tar.gz &&\
    tar xf v$SUNDIALS_VERSION.tar.gz && rm -f v$SUNDIALS_VERSION.tar.gz && \
    cd sundials-$SUNDIALS_VERSION && \
    mkdir build && cd build && \
    cmake \
        -DCMAKE_INSTALL_PREFIX=$HOME/libs/sundials-$SUNDIALS_VERSION \
        -DLAPACK_ENABLE=ON \
        -DMPI_ENABLE=ON .. && \
    make -j4 && make install && \
    cd $HOME && \
    rm -rf sundials-$SUNDIALS_VERSION
ENV SUNDIALS_DIR $HOME/libs/sundials-$SUNDIALS_VERSION

ARG D2K_START=3c27de2

RUN \
    wget https://github.com/mathLab/deal2lkit/archive/master.zip && \
    unzip master.zip && rm master.zip && \
    cd deal2lkit-master && mkdir build && cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=$HOME/libs/deal2lkit-master \
          -DCMAKE_BUILD_TYPE=Debug \
	  -DD2K_COMPONENT_DOCUMENTATION=OFF \
	  -GNinja \
          ../ && \
    ninja -j4 && ninja install && \
    cd $HOME && rm -rf deal2lkit-master

ENV D2K_DIR $HOME/libs/deal2lkit-master
ENV DEAL2LKIT_DIR $HOME/libs/deal2lkit-master
