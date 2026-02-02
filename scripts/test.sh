#!/bin/bash

qsup scripts/test_seq_omp.sh
qsup scripts/test_omp2.sh
qsub scripts/test_gapbs.sh