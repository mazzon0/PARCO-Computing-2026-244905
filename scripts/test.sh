#!/bin/bash

qsub scripts/test_seq_omp.pbs
qsub scripts/test_omp2.pbs
qsub scripts/test_mpi.pbs
qsub scripts/test_gapbs.pbs