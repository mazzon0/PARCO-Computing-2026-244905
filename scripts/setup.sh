#!/bin/bash

set -euo pipefail

# Load GCC and G++
echo "Loading GCC and MPI..."
module load gcc91
gcc () { gcc -9.1.0 "$@"; }
module load mpich-3.2.1-gcc-9.1.0
echo "Loaded GCC"
echo ""

# Compile
echo "Compiling..."
make clean
make
echo "Compilation finished"
echo ""

# Download graphs
echo "Downloading graphs..."
scripts/download.sh stanford
scripts/download.sh google
echo "Graphs downloaded"
echo ""

# Unload modules
module unload gcc91
module unload module load mpich-3.2.1-gcc-9.1.0

# TODO add other implementations to test