#include <mpi.h>
#include <stdbool.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "common/csr_utils.h"

double *e;
double *vec0;
double *vec1;
uint64_t *dangling;
double *global_vec;
int *counts;
int *displs;

bool pagerank(const csr_matrix_t *const mat, double **rankv, int rank, int n_proc);
void cleanup(csr_matrix_t *const mat);
int load_csr_mpi(const char *filename, csr_matrix_t *mat, int rank, int n);
void find_dangling(const csr_matrix_t *const mat, uint64_t *const out, uint64_t *out_size, int rank, int n_proc);
void matvec_mul(const csr_matrix_t *const mat, const double *const vec, double *const out, int n_proc);
double rank_loss(const uint64_t *const dangling_indices, const uint64_t dangling_size, const double *const rank);
double random_surfer(double *new_rank, double *last_rank, double *e, uint64_t local_size, uint64_t global_size, double random_jump_prob, double teleport);

int main(int argc, char **argv) {
    // Check errors
    if (argc < 2 || argc > 3) {
        fprintf(stderr, "Usage: %s path/to/web_graph.csr [tries]\n", argv[0]);
        return -1;
    }

    MPI_Init(&argc, &argv);
    int rank, n_proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);

    // Parse tries (default to 10 if not provided)
    int tries = 10; 
    if (argc == 3) {
        tries = atoi(argv[2]);
        if (tries <= 0) {
            fprintf(stderr, "Error: Number of tries must be a positive integer.\n");
            MPI_Finalize();
            return -1;
        }
    }
    
    // Load web graph
    csr_matrix_t web;
    if (load_csr_mpi(argv[1], &web, rank, n_proc) != 0) {
        MPI_Finalize();
        return -1;
    }

    // Preallocate buffers
    uint64_t local_size = web.n_rows;
    uint64_t global_size = web.n_columns;

    e = malloc(sizeof(double) * local_size);
    vec0 = malloc(sizeof(double) * local_size);
    vec1 = malloc(sizeof(double) * local_size);
    dangling = malloc(sizeof(uint64_t) * local_size);
    global_vec = malloc(sizeof(double) * global_size);
    counts = malloc(sizeof(int) * n_proc);
    displs = malloc(sizeof(int) * n_proc);

    int local_n_int = (int)local_size;
    MPI_Allgather(&local_n_int, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);
    displs[0] = 0;
    for (int i = 1; i < n_proc; i++) {
        displs[i] = displs[i - 1] + counts[i - 1];
    }

    // Run PageRank
    double *rankv = malloc(sizeof(double) * local_size);
    for (int i = 0; i < tries; i++) {
        pagerank(&web, &rankv, rank, n_proc);
    }

    free(rankv);
    free(e);
    free(vec0);
    free(vec1);
    free(dangling);
    free(global_vec);
    free(counts);
    free(displs);
    cleanup(&web); 

    MPI_Finalize();
    return 0;
}

bool pagerank(const csr_matrix_t *const mat, double **rankv, int rank, int n_proc) {
    const double EPSILON = 1e-6;
    const double GLOBAL_RANK = 1.0;
    const uint32_t MAX_ITERATIONS = 1000;
    const double RANDOM_JUMP_PROB = 0.15;
    
    uint64_t local_size = mat->n_rows;
    uint64_t global_size = mat->n_columns;

    // Initial distribution: uniform
    double init_val = GLOBAL_RANK / global_size;
    for (uint64_t i = 0; i < local_size; i++) {
        e[i] = init_val;
    }

    // Initialization
    double *last_rank = vec0;
    double *new_rank = vec1;
    memcpy(vec0, e, sizeof(double) * local_size);
    uint64_t dangling_size;
    find_dangling(mat, dangling, &dangling_size, rank, n_proc);

    double delta;
    uint32_t iteration = 0;

    do {
        // Matrix-Vector multiplication
        matvec_mul(mat, last_rank, new_rank, n_proc);
        
        // Compute lost rank (dangling pages)
        double teleport = rank_loss(dangling, dangling_size, new_rank);
        
        // Random surfer + compute delta
        delta = random_surfer(new_rank, last_rank, e, local_size, global_size, RANDOM_JUMP_PROB, teleport);

        // Update data of last iteration
        double *aux = new_rank;
        new_rank = last_rank;
        last_rank = aux;

        iteration++;
        // Check convergence
    } while(delta > EPSILON && iteration < MAX_ITERATIONS);

    // Store result
    memcpy(*rankv, last_rank, sizeof(double) * local_size);

    return true;
}

void cleanup(csr_matrix_t *const mat) {
    free(mat->values);
    free(mat->columns);
    free(mat->row_ptrs);
}

int load_csr_mpi(const char *filename, csr_matrix_t *mat, int rank, int n) {
    MPI_File fh;
    if (MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh)) {
        return -1;
    }    

    // Read header
    uint64_t header[3]; // n_rows, n_columns, nnz
    MPI_File_read_at_all(fh, 0, header, 3, MPI_UINT64_T, MPI_STATUS_IGNORE);
    mat->n_rows = header[0];
    mat->n_columns = header[1];
    mat->nnz = header[2];

    // Read row_ptrs
    uint64_t local_n = mat->n_rows / n;
    uint64_t start_row = rank * local_n;
    if (rank == n - 1) local_n = mat->n_rows - start_row;
    
    mat->row_ptrs = malloc((local_n + 1) * sizeof(uint64_t));
    MPI_Offset row_ptr_offset = 24 + (mat->nnz * sizeof(double)) + (mat->nnz * sizeof(uint64_t)) + (start_row * sizeof(uint64_t));
    MPI_File_read_at_all(fh, row_ptr_offset, mat->row_ptrs, local_n, MPI_UINT64_T, MPI_STATUS_IGNORE);

    uint64_t next_ptr;
    if (start_row + local_n < mat->n_rows) {
        MPI_File_read_at(fh, row_ptr_offset + (local_n * sizeof(uint64_t)), &next_ptr, 1, MPI_UINT64_T, MPI_STATUS_IGNORE);
    } else {
        next_ptr = mat->nnz;
    }
    mat->row_ptrs[local_n] = next_ptr;

    // Find local NNZ
    uint64_t local_nnz = mat->row_ptrs[local_n] - mat->row_ptrs[0];
    uint64_t nnz_offset_count = mat->row_ptrs[0];

    mat->values = malloc(local_nnz * sizeof(double));
    mat->columns = malloc(local_nnz * sizeof(uint64_t));

    MPI_File_read_at_all(fh, 24 + (nnz_offset_count * sizeof(double)), mat->values, local_nnz, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_read_at_all(fh, 24 + (mat->nnz * sizeof(double)) + (nnz_offset_count * sizeof(uint64_t)), mat->columns, local_nnz, MPI_UINT64_T, MPI_STATUS_IGNORE);

    uint64_t base_offset = mat->row_ptrs[0];
    for(uint64_t i = 0; i <= local_n; i++) mat->row_ptrs[i] -= base_offset;
    mat->n_rows = local_n;
    mat->nnz = local_nnz;

    MPI_File_close(&fh);
    return 0;
}

void find_dangling(const csr_matrix_t *const mat, uint64_t *const out, uint64_t *out_size, int rank, int n_proc) {
    uint64_t global_size = mat->n_columns;
    
    uint32_t *local_counts = calloc(global_size, sizeof(uint32_t));
    uint32_t *global_counts = malloc(global_size * sizeof(uint32_t));

    // Count local links per column
    for (uint64_t i = 0; i < mat->nnz; i++) {
        local_counts[mat->columns[i]]++;
    }

    MPI_Allreduce(local_counts, global_counts, (int)global_size, MPI_UINT32_T, MPI_SUM, MPI_COMM_WORLD);

    *out_size = 0;
    uint64_t start_row = displs[rank];
    
    // Construct the array with local dangling pages
    for (uint64_t i = 0; i < mat->n_rows; i++) {
        if (global_counts[start_row + i] == 0) {
            out[*out_size] = i;
            (*out_size)++;
        }
    }

    free(local_counts);
    free(global_counts);
}

void matvec_mul(const csr_matrix_t *const mat, const double *const vec, double *const out, int n_proc) {
    MPI_Allgatherv(vec, (int)mat->n_rows, MPI_DOUBLE, global_vec, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

    for (uint64_t j = 0; j < mat->n_rows; j++) {
        double sum = 0.0;
        for (uint64_t index = mat->row_ptrs[j]; index < mat->row_ptrs[j + 1]; index++) {
            sum += mat->values[index] * global_vec[mat->columns[index]];
        }
        out[j] = sum;
    }
}

double rank_loss(const uint64_t *const dangling_indices, const uint64_t dangling_size, const double *const rank) {
    double local_sum = 0.0;
    for (uint64_t i = 0; i < dangling_size; i++) {
        local_sum += rank[dangling_indices[i]];
    }
    
    double global_sum = 0.0;
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_sum;
}

double random_surfer(double *new_rank, double *last_rank, double *e, uint64_t local_size, uint64_t global_size, double random_jump_prob, double teleport) {
    double local_delta = 0.0;

    for (uint64_t i = 0; i < local_size; i++) {
        // Random Surfer update
        double v = (1.0 - random_jump_prob) * new_rank[i] + random_jump_prob * e[i];
        v += teleport / global_size;
        new_rank[i] = v;
        
        // Compute delta
        local_delta += fabs(v - last_rank[i]);
    }
    
    double global_delta = 0.0;
    MPI_Allreduce(&local_delta, &global_delta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_delta;
}
