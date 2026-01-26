#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s web-Stanford.txt output.csr\n", argv[0]);
        return 1;
    }

    const char *input = argv[1];
    const char *output = argv[2];

    FILE *f = fopen(input, "r");
    if (!f) {
        perror("fopen input");
        return 1;
    }

    // Count edges and maximum node
    uint64_t from, to;
    uint64_t max_node = 0;
    uint64_t nnz = 0;
    char line[256];

    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#') continue;
        if (sscanf(line, "%lu %lu", &from, &to) == 2) {
            if (from > max_node) max_node = from;
            if (to > max_node) max_node = to;
            nnz++;
        }
    }

    uint64_t n_rows = max_node;
    uint64_t n_columns = max_node;

    // Count per-row and per-column links
    uint64_t *row_counts = calloc(n_rows, sizeof(uint64_t));
    uint64_t *col_elements = calloc(n_columns, sizeof(uint64_t));

    rewind(f);
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#') continue;
        if (sscanf(line, "%lu %lu", &from, &to) == 2) {
            row_counts[from]++;       // outgoing edges
            col_elements[to]++;       // incoming edges
        }
    }

    // Build row_ptrs
    uint64_t *row_ptrs = calloc(n_rows, sizeof(uint64_t));
    uint64_t sum = 0;
    for (uint64_t i = 0; i < n_rows; i++) {
        row_ptrs[i] = sum;
        sum += row_counts[i];
    }

    // Fill CSR arrays
    double   *values  = malloc(nnz * sizeof(double));
    uint64_t *columns = malloc(nnz * sizeof(uint64_t));
    uint64_t *offsets = calloc(n_rows, sizeof(uint64_t));

    rewind(f);
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#') continue;
        if (sscanf(line, "%lu %lu", &from, &to) == 2) {
            uint64_t row = from;
            uint64_t idx = row_ptrs[row] + offsets[row];

            // Safety check
            if (idx >= nnz) {
                fprintf(stderr, "ERROR: idx %lu >= nnz %lu (from=%lu, row=%lu, offsets[row]=%lu)\n",
                        idx, nnz, from, row, offsets[row]);
                exit(1);
            }

            columns[idx] = to;
            values[idx]  = 1.0;
            offsets[row]++;
        }
    }

    free(offsets);

    // Normalize values
    for (uint64_t i = 0; i < nnz; i++) {
        values[i] /= col_elements[columns[i]];
    }

    // Write binary
    FILE *out = fopen(output, "wb");
    if (!out) {
        perror("fopen output");
        return 1;
    }

    fwrite(&n_rows,    sizeof(uint64_t), 1, out);
    fwrite(&n_columns, sizeof(uint64_t), 1, out);
    fwrite(&nnz,       sizeof(uint64_t), 1, out);
    fwrite(values,     sizeof(double),   nnz, out);
    fwrite(columns,    sizeof(uint64_t), nnz, out);
    fwrite(row_ptrs,   sizeof(uint64_t), n_rows, out);

    fclose(out);

    free(values);
    free(columns);
    free(row_ptrs);

    printf("CSR written successfully:\n");
    printf("\trows    = %lu\n", n_rows);
    printf("\tcolumns = %lu\n", n_columns);
    printf("\tnnz     = %lu\n", nnz);

    return 0;
}
