docker run --privileged -t --rm -P -v `pwd`:/home/dealii/deal2lkit:rw mathlab/deal2lkit-base:v8.5.0 /bin/sh -c "$@"
