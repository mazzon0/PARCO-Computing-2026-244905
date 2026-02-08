#include "csr_utils.h"

int load_csr(const char *filename, csr_matrix_t *mat) {
    FILE *f = fopen(filename, "rb");
    if (!f) {
        perror("fopen");
        return -1;
    }

    // Read n_rows, n_columns, nnz
    if (fread(&mat->n_rows, sizeof(uint64_t), 1, f) != 1 ||
        fread(&mat->n_columns, sizeof(uint64_t), 1, f) != 1 ||
        fread(&mat->nnz, sizeof(uint64_t), 1, f) != 1) {
        perror("fread metadata");
        fclose(f);
        return -1;
    }

    // Allocate memory
    mat->values = malloc(mat->nnz * sizeof(float));
    mat->columns = malloc(mat->nnz * sizeof(uint64_t));
    mat->row_ptrs = malloc((mat->n_rows + 1) * sizeof(uint64_t));

    if (!mat->values || !mat->columns || !mat->row_ptrs) {
        perror("malloc");
        fclose(f);
        free(mat->values);
        free(mat->columns);
        free(mat->row_ptrs);
        return -1;
    }

    // Read values
    if (fread(mat->values, sizeof(float), mat->nnz, f) != mat->nnz) {
        perror("fread values");
        fclose(f);
        free(mat->values); free(mat->columns); free(mat->row_ptrs);
        return -1;
    }

    // Read columns
    if (fread(mat->columns, sizeof(uint64_t), mat->nnz, f) != mat->nnz) {
        perror("fread columns");
        fclose(f);
        free(mat->values); free(mat->columns); free(mat->row_ptrs);
        return -1;
    }

    // Read row_ptrs
    if (fread(mat->row_ptrs, sizeof(uint64_t), mat->n_rows + 1, f) != mat->n_rows + 1) {
        perror("fread row_ptrs");
        fclose(f);
        free(mat->values); free(mat->columns); free(mat->row_ptrs);
        return -1;
    }
    mat->row_ptrs[mat->n_rows] = mat->nnz;

    fclose(f);
    return 0;
}
