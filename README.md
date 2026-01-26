# Parallel PageRank

This repository implements the PageRank algorithm and tries multiple methods to parallelize it. This is a project for the "Introduction to Parallel Computing" course at University of Trento.

## Working on a remote computer
Coming soon.

## Working on a local computer

### Installation
This project requires a Linux environment with some basic development packages: `gcc`, `gzip`, `bash`, `wget`.

On Ubuntu and Debian you can run this command to install the necessary packages.
```bash
sudo apt-get update && sudo apt-get install build-essential gzip wget
```

Download the repo
```bash
git clone git@github.com:mazzon0/PARCO-Computing-2026-244905.git
```

### Compilation
The project is composed of 4 programs:
- `seq` is a sequential version of PageRank
- `omp` is a parallel version of PageRank (using OpenMP)
- `mpi` is a parallel version of PageRank (using MPI)
- `conv` is a program necessary to convert downloaded dataset to CSR

You can compile everything with `make` or compile just the needed programs.
```bash
make seq
make omp
make mpi
make conv
```

The executables will be places in `bin/`.

### Datasets
Some commonly used datasets can be downloaded and converted automatically with a script.
```bash
./download.sh dataset_name
```

Run ```./download.sh``` to see which datasets are available. Datasets will be downloaded into `datasets/`.

### Execution
Execute the sequantial version with
```bash
bin/seq datasets/dataset_name.csr
```

OMP and MPI are coming soon.
