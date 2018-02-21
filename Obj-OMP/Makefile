# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -fopenmp
LDLIBS=-lm

# list of source files
SRC=ljmd-omp.c
VPATH=../src
OBJ=$(SRC:%.c=%.o)

default: ../ljmd-omp.x

clean:
	rm -f *.mod *.o ../ljmd-omp.x

# linker rule
../ljmd-omp.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

# dependencies
ljmd-omp.o: ljmd-omp.c
