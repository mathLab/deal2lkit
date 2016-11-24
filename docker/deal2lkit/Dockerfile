FROM mathlab/deal2lkit-base

MAINTAINER luca.heltai@gmail.com

ARG BUILD_TYPE=DebugRelease
ARG VER=master
ARG KEEP_SOURCE=true

RUN \
    git clone https://github.com/mathLab/deal2lkit/ deal2lkit-$VER-src && \
    cd deal2lkit-$VER-src && git checkout $VER && \
    mkdir build && cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=$HOME/libs/deal2lkit-$VER\
          -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
	  -DD2K_COMPONENT_DOCUMENTATION=OFF \
	  -GNinja \
          ../ && \
    ninja -j4 && ninja install && \
    cd $HOME/deal2lkit-$VER-src/ && \ 
    ($KEEP_SOURCE && rm -rf build .git) || \
    cd $HOME && rm -rf deal2lkit-$VER-src

ENV D2K_DIR $HOME/libs/deal2lkit-master
ENV DEAL2LKIT_DIR $HOME/libs/deal2lkit-master
