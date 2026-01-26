CFLAGS = -O2 -Wall

all: seq omp mpi conv

seq:
	mkdir -p bin
	gcc $(CFLAGS) src/pagerank_seq.c -Iinclude -o bin/seq

omp:
	mkdir -p bin
	gcc -fopenmp $(CFLAGS) src/pagerank_omp.c -Iinclude -o bin/omp

mpi:
	mkdir -p bin
	mpicc $(CFLAGS) src/pagerank_mpi.c -Iinclude -o bin/mpi

conv:
	mkdir -p bin
	gcc $(CFLAGS) src/converter.c -Iinclude -o bin/converter

clean:
	rm -f bin/*
