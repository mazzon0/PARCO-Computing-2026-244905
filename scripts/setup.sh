#!/bin/bash

set -euo pipefail

# Load GCC and G++
module load gcc91
module load g++91
gcc () { gcc -9.1.0 "$@"; }
g++ () { g++ -9.1.0 "$@"; }
gcc --version
g++ --version

# Compile
make clean
make

# Download graphs
scripts/download.sh stanford
scripts/download.sh google

# TODO add other implementations to test