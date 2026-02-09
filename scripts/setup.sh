#!/bin/bash

set -euo pipefail

# Load GCC and G++
echo "Loading GCC and MPI..."
module load gcc91
gcc () { gcc -9.1.0 "$@"; }
module load mpich-3.2.1--gcc-9.1.0
echo "Loaded GCC"
echo ""

# Setup Python environment
scripts/setup_python.sh

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

# Synthetic graphs
echo "Generating synthetic graphs..."
scripts/.venv/bin/python scripts/gen.py -n 100000 -z 79057 -f datasets/100k
scripts/.venv/bin/python scripts/gen.py -n 200000 -z 111804 -f datasets/200k
scripts/.venv/bin/python scripts/gen.py -n 400000 -z 158114 -f datasets/400k
scripts/.venv/bin/python scripts/gen.py -n 800000 -z 223607 -f datasets/800k
scripts/.venv/bin/python scripts/gen.py -n 1600000 -z 316228 -f datasets/1600k
scripts/.venv/bin/python scripts/gen.py -n 3200000 -z 447214 -f datasets/3200k
scripts/.venv/bin/python scripts/gen.py -n 6400000 -z 632456 -f datasets/6400k
echo "Graphs generated."
echo ""

# Unload modules
module unload gcc91
module unload module load mpich-3.2.1--gcc-9.1.0

# TODO add other implementations to test