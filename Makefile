CFLAGS = -std=c99 -O3 -Wall -march=native
COMMON_SRC = src/common/csr_utils.c
LDFLAGS = -lm -lrt

all: seq omp omp2 mpi conv

seq:
	mkdir -p bin
	gcc $(CFLAGS) src/pagerank_seq.c $(COMMON_SRC) -Iinclude -o bin/seq $(LDFLAGS)

omp:
	mkdir -p bin
	gcc -fopenmp $(CFLAGS) src/pagerank_omp.c $(COMMON_SRC) -Iinclude -o bin/omp $(LDFLAGS)

omp2:
	mkdir -p bin
	gcc -fopenmp $(CFLAGS) src/pagerank_omp.c $(COMMON_SRC) -Iinclude -o bin/omp2 -DMERGE_VECTOR_OPERATIONS $(LDFLAGS)

mpi:
	mkdir -p bin
	mpicc $(CFLAGS) src/pagerank_mpi.c -Iinclude -o bin/mpi $(LDFLAGS)

conv:
	mkdir -p bin
	gcc $(CFLAGS) src/converter.c -Iinclude -o bin/converter

clean:
	rm -rf bin/*
