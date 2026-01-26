CFLAGS = -O2 -Wall

all: seq omp mpi

seq:
	mkdir -p bin
	gcc $(CFLAGS) src/pagerank_seq.c -Iinclude -o bin/seq

omp:
	gcc -fopenmp $(CFLAGS) src/pagerank_omp.c -Iinclude -o bin/omp

mpi:
	mpicc $(CFLAGS) src/pagerank_mpi.c -Iinclude -o bin/mpi

clean:
	rm -f bin/*
