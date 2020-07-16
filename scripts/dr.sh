docker run -t --rm -P -v `pwd`:/home/dealii/deal2lkit:rw dealii/dealii:v9.2.0-focal /bin/sh -c "cd deal2lkit; $@"
