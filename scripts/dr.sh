docker run -i -t --rm -P -v `pwd`:/home/dealii/deal2lkit:rw dealii/dealii:v9.1.1-gcc-mpi-fulldepsspack-debugrelease /bin/sh -c "cd deal2lkit; $@"
