#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include "common/csr_utils.h"

bool pagerank(const csr_matrix_t *const mat, double **rank);
bool pagerank_original(const csr_matrix_t *const mat, double **rank);
void cleanup(csr_matrix_t *mat, double *rank);

void matvec_mul(const csr_matrix_t *const mat, const double *const vec, double *const out);
void vec_add_scalar(const double *const vec1, double scalar, double *const out, uint64_t size);
void vec_diff(const double *const vec1, const double *const vec2, double *const out, uint64_t size);
void linear_comb(const double *const vec1, const double *const vec2, const double a1, const double a2, double *const out, uint64_t size);
double l1_norm(const double *const vec, uint64_t size);

void find_dangling(const csr_matrix_t *const mat, uint64_t *const out, uint64_t *out_size);
double rank_loss(const uint64_t *const dangling_indices, const uint64_t dangling_size, const double *const rank, const uint64_t size);


int main(int argc, char **argv) {
    // Check errors
    if (argc != 2) {
        fprintf(stderr, "Usage: %s path/to/web_graph.csr\n", argv[0]);
        return -1;
    }

    // Load web graph
    csr_matrix_t web;
    if (load_csr(argv[1], &web) != 0) {
        return -1;
    }

    // Run PageRank
    double *rank = malloc(sizeof(double) * web.n_rows);
    if (!pagerank(&web, &rank)) {
        cleanup(&web, rank);
        return -1;
    }

    // Cleanup
    cleanup(&web, rank);

    return 0;
}


/**
 * PageRank implementation
 */
bool pagerank(const csr_matrix_t *const mat, double **rank) {
    if (mat->n_rows != mat->n_columns) {
        fprintf(stderr, "ERROR: n_rows (%" PRIu64 ") must be equal to n_columns (%" PRIu64 ")", mat->n_rows, mat->n_columns);
        return false;
    }

    const double EPSILON = 1e-7;
    const double GLOBAL_RANK = 1.0;
    const uint32_t MAX_ITERATIONS = 1000;
    const double RANDOM_JUMP_PROB = 0.15;
    
    uint64_t size = mat->n_rows;

    // Initial distribution: uniform
    double *e = malloc(sizeof(double) * size);
    for (uint64_t i = 0; i < size; i++) {
        e[i] = GLOBAL_RANK / size;
    }

    // Resource allocation
    double *vec0 = malloc(sizeof(double) * size);
    double *vec1 = malloc(sizeof(double) * size);
    double *diff = malloc(sizeof(double) * size);
    uint64_t *dangling = malloc(sizeof(uint64_t) * size);

    // Initialization
    double *last_rank = vec0;
    double *new_rank = vec1;
    memcpy(vec0, e, sizeof(double) * size);
    uint64_t dangling_size;
    find_dangling(mat, dangling, &dangling_size);

    double delta;
    uint32_t iteration = 0;

    do {
        // Matrix-Vector multiplication
        matvec_mul(mat, last_rank, new_rank);
        
        // Compute lost rank (dangling pages)
        double teleport = rank_loss(dangling, dangling_size, new_rank, size);
        
#       ifndef MERGE_VECTOR_OPERATIONS
        
        // Random surfer
        linear_comb(new_rank, e, 1.0 - RANDOM_JUMP_PROB, RANDOM_JUMP_PROB, new_rank, size);
        vec_add_scalar(new_rank, teleport / (double)size, new_rank, size);

        // Compute delta
        vec_diff(new_rank, last_rank, diff, size);
        delta = l1_norm(diff, size);

#       else

        delta = 0.0;

        #pragma omp parallel for reduction(+:delta)
        for (uint64_t i = 0; i < size; i++) {
            double v = (1.0 - RANDOM_JUMP_PROB) * new_rank[i] + RANDOM_JUMP_PROB * e[i];
            v += teleport / size;
            diff[i] = v - last_rank[i];
            delta += fabs(v - last_rank[i]);
            new_rank[i] = v;
        }

#       endif

        // Update data of last iteration
        double *aux = new_rank;
        new_rank = last_rank;
        last_rank = aux;

        iteration++;

        // Check convergence
    } while(delta > EPSILON && iteration < MAX_ITERATIONS);

    // Store result
    memcpy(*rank, last_rank, sizeof(double) * size);

    // Cleanup
    free(e);
    free(vec0);
    free(vec1);
    free(diff);
    free(dangling);

    return true;
}

/**
 * Frees the resources used by the main function
 */
void cleanup(csr_matrix_t *mat, double *rank) {
    free(rank);
    free(mat->values);
    free(mat->columns);
    free(mat->row_ptrs);
}


// Helper functions

void matvec_mul(const csr_matrix_t *const mat, const double *const vec, double *const out) {
    #pragma omp parallel for schedule(static, 512)
    for (uint64_t j = 0; j < mat->n_rows; j++) {
        double sum = 0.0;

        for (uint64_t index = mat->row_ptrs[j]; index < mat->row_ptrs[j + 1]; index++) {
            sum += mat->values[index] * vec[mat->columns[index]];
        }

        out[j] = sum;
    }
}

void vec_add_scalar(const double *const vec1, double scalar, double *const out, uint64_t size) {
    #pragma omp parallel for
    for (uint64_t i = 0; i < size; i++) {
        out[i] = vec1[i] + scalar;
    }
}

void vec_diff(const double *const vec1, const double *const vec2, double *const out, uint64_t size) {
    #pragma omp parallel for
    for (uint64_t i = 0; i < size; i++) {
        out[i] = vec1[i] - vec2[i];
    }
}

void linear_comb(const double *const vec1, const double *const vec2, const double a1, const double a2, double *const out, uint64_t size) {
    #pragma omp parallel for
    for (uint64_t i = 0; i < size; i++) {
        out[i] = vec1[i] * a1 + vec2[i] * a2;
    }
}

double l1_norm(const double *const vec, uint64_t size) {
    double sum = 0.0;

    #pragma omp parallel for reduction(+:sum)
    for (uint64_t i = 0; i < size; i++) {
        sum += fabs(vec[i]);
    }

    return sum;
}

void find_dangling(const csr_matrix_t *const mat, uint64_t *const out, uint64_t *out_size) {
    // Count links of each page (using temporarily the out buffer)
    for (uint64_t i = 0; i < mat->nnz; i++) {
        out[mat->columns[i]] += 1;
    }

    // Set the ids of the dangling pages
    *out_size = 0;
    for (uint64_t i = 0; i < mat->n_columns; i++) {
        if (out[i] == 0) {
            out[*out_size] = i;
            (*out_size)++;
        }
    }
}

double rank_loss(const uint64_t *const dangling_indices, const uint64_t dangling_size, const double *const rank, const uint64_t size) {
    double sum = 0.0;

    #pragma omp parallel for reduction(+:sum)
    for (uint64_t i = 0; i < dangling_size; i++) {
        sum += rank[dangling_indices[i]];
    }

    return sum;
}
