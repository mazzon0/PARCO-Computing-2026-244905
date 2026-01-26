#!/bin/bash

dataset=$1

case $dataset in
    "stanford")
        mkdir -p datasets/
        cd datasets/

        echo "Downloading dataset..."
        wget https://snap.stanford.edu/data/web-Stanford.txt.gz

        echo "Decompressing dataset..."
        gunzip web-Stanford.txt.gz
        cd ..

        echo "Converting to CSR"
        if [ ! -x bin/converter ]; then
            echo "Error: bin/converter not found or not executable."
            echo "Run 'make conv' and try again."
            exit 1
        fi
        bin/converter datasets/web-Stanford.txt datasets/stanford.csr
        rm datasets/web-Stanford.txt

        echo "Done. Run PageRank with 'bin/seq datasets/stanford.csr'"
        ;;

    *)
        echo "Usage: ./download.sh {stanford}"
        ;;
esac
