# Parallel PageRank

This repository implements the PageRank algorithm and tries multiple methods to parallelize it.
This is a project for the "Introduction to Parallel Computing" course at University of Trento.

PageRank, developed by Google founders Sergey Brin and Larry Page, ranks web pages based on their link structure to measure their importance. 
The algorithm involves iterative computations over sparse graphs represented by large matrices. 
Parallelizing PageRank can dramatically reduce execution time by distributing computations across multiple cores or nodes, 
making it an interesting problem to explore parallel computing frameworks such as OpenMP and MPI.

## Working on a remote cluster
This section includes all the information to work with this repo in a remote cluster with the PBS scheduler.

### Quick Start (UniTN cluster)
These section contains the instructions to quickly setup the project on the cluster `hpc.unitn.it`.
To reproduce the project in other environments, please check out the following sections.

First, start an interactive session on a compute node.
```bash
qsub -I -q short_cpuQ
```

Download the repository.
```bash
git clone git@github.com:mazzon0/PARCO-Computing-2026-244905.git
cd PARCO-Computing-2026-244905
```

Compile the source code and download the datasets.
```bash
scripts/setup.sh
```

Exit the interactive session and go back to a head node.
```bash
exit
```

Now the setup has completed, and you are ready to schedule the PBS job.
```bash
cd PARCO-Computing-2026-244905
qsub scripts/test.pbs
```

### TODO

## Working on a local computer

### Installation
This project requires a Linux environment with some basic development packages: `gcc`, `make`, `gzip`, `bash`, `wget`.

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
- `omp` is a parallel version of PageRank with shared memory (using OpenMP)
- `mpi` is a parallel version of PageRank with distributed memory (using MPI)
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
Some commonly used datasets can be downloaded and converted automatically in the `datasets/` directory.
```bash
scripts/download.sh dataset_name
```

Run `scripts/download.sh` to see which datasets are available.

### Execution
Execute the sequential version with
```bash
bin/seq datasets/dataset_name.csr
```

Execute the OpenMP version with
```bash
bin/omp datasets/dataset_name.csr
```

MPI is coming soon.
