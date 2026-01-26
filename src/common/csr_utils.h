#ifndef CSR_UTILS
#define CSR_UTILS

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    uint64_t n_rows;
    uint64_t n_columns;
    uint64_t nnz;
    double *values;
    uint64_t *columns;
    uint64_t *row_ptrs;
} csr_matrix_t;

int load_csr(const char *filename, csr_matrix_t *mat);

#endif
