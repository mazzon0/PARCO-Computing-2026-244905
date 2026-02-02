#!/bin/bash
set -e

mkdir -p comparison
cd comparison
echo "Downloading GAPBS..."
git clone https://github.com/sbeamer/gapbs.git
echo "Download finished"
echo ""

cd gapbs
echo "Compiling PageRank (SpMV)..."
make pr_spmv
echo "Compilation finished"