#!/bin/bash

program=$1
times=$2

case $program in
    "seq")
        ;;
    
    "omp")
        export OMP_NUM_THREADS=$3
        ;;

    *)
        echo "Usage 'scripts/test.sh {seq|omp} {num_iterations}'"
        exit 1
esac

declare -a datasets=("stanford" "google")

for ((i=0; i<$times; i++)); do
    echo "iteration: $i"
    for dataset in "${datasets[@]}"; do
        perf stat -e cycles,instructions,cache-references,cache-misses bin/$program datasets/$dataset.csr 2>> perf.txt
    done
done

echo "Done."
