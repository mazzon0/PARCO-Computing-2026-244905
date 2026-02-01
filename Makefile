CFLAGS = -std=c99 -O2 -Wall
COMMON_SRC = src/common/csr_utils.c

all: seq omp conv

seq:
	mkdir -p bin
	gcc $(CFLAGS) src/pagerank_seq.c $(COMMON_SRC) -Iinclude -o bin/seq

omp:
	mkdir -p bin
	gcc -fopenmp $(CFLAGS) src/pagerank_omp.c $(COMMON_SRC) -Iinclude -o bin/omp

mpi:
	mkdir -p bin
	mpicc $(CFLAGS) src/pagerank_mpi.c $(COMMON_SRC) -Iinclude -o bin/mpi

conv:
	mkdir -p bin
	gcc $(CFLAGS) src/converter.c -Iinclude -o bin/converter

clean:
	rm -f bin/*
