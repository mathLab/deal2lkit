FROM heltai/dealii-release

MAINTAINER luca.heltai@gmail.com

# deal2lkit master image
ENV D2K_START 3c27de2
RUN \
    wget https://github.com/mathLab/deal2lkit/archive/master.zip && \
    unzip master.zip && rm master.zip && \
    cd deal2lkit-master && mkdir build && cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=$HOME/libs/deal2lkit-master \
          -DCMAKE_BUILD_TYPE=Release \
	  -DD2K_COMPONENT_DOCUMENTATION=OFF \
	  -GNinja \
          ../ && \
    ninja -j4 && ninja install && \
    cd $HOME && rm -rf deal2lkit-master

ENV D2K_DIR $HOME/libs/deal2lkit-master
ENV DEAL2LKIT_DIR $HOME/libs/deal2lkit-master

# Need this for Travis!
USER root
