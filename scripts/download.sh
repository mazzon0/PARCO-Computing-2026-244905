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
        grep -v '^#' datasets/web-Stanford.txt | awk '{print $1" "$2}' > datasets/stanford.el
        rm datasets/web-Stanford.txt

        echo "Done. Run PageRank with 'bin/seq datasets/stanford.csr'"
        ;;

    "google")
        mkdir -p datasets/
        cd datasets/

        echo "Downloading dataset..."
        wget https://snap.stanford.edu/data/web-Google.txt.gz

        echo "Decompressing dataset..."
        gunzip web-Google.txt.gz
        cd ..

        echo "Converting to CSR"
        if [ ! -x bin/converter ]; then
            echo "Error: bin/converter not found or not executable."
            echo "Run 'make conv' and try again."
            exit 1
        fi
        bin/converter datasets/web-Google.txt datasets/google.csr
        grep -v '^#' datasets/web-Google.txt | awk '{print $1" "$2}' > datasets/google.el
        rm datasets/web-Google.txt

        echo "Done. Run PageRank with 'bin/seq datasets/google.csr'"
        ;;

    *)
        echo "Usage: scripts/download.sh {google|stanford}"
        echo "  google:     5,105,039 links"
        echo "  stanford:   2,312,497 links"
        ;;
esac
