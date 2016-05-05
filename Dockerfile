FROM heltai/dealii

MAINTAINER luca.heltai@gmail.com

# deal2lkit repo
ENV D2K_START 9762009
RUN git clone https://github.com/mathLab/deal2lkit.git

ENV D2K_START 9762009
#build
RUN cd deal2lkit && \
    mkdir build && cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=$HOME/deal2lkit-inst \
          -DCMAKE_BUILD_TYPE=Debug \
	  -DD2K_COMPONENT_DOCUMENTATION=OFF \
	  -GNinja \
          ../ && \
    ninja -j4 && ninja install && \
    cd .. && rm -rf build

ENV D2K_DIR $HOME/deal2lkit-inst
ENV DEAL2LKIT_DIR $HOME/deal2lkit-inst
