#! /bin/sh
cmd() {
  echo "+ $@"
  eval "$@"
}

CC=gcc-13
CFLAGS=-fopenmp
# echo $CC $CFLAGS $1
cmd $CC $CFLAGS $1
for i in 1 2 4 8
    do cmd OMP_NUM_THREADS=$i ./a.out
done
