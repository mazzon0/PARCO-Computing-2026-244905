#define _POSIX_C_SOURCE 199309L
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include "common/csr_utils.h"
#include "common/time_utils.h"

typedef struct {
    uint64_t original_id;
    uint32_t degree;
} node_degree_t;

double *execution_times;

float *e;
float *vec0;
float *vec1;
float *diff;
uint64_t *dangling;

bool pagerank(const csr_matrix_t *const mat, float **rank);
bool pagerank_original(const csr_matrix_t *const mat, float **rank);
void cleanup(csr_matrix_t *mat, float *rank);
int compare_doubles(const void *a, const void *b);
double get_average(double *times, int n);
double get_median(double *times, int n);

void matvec_mul(const csr_matrix_t *const mat, const float *const vec, float *const out);
void vec_add_scalar(const float *const vec1, float scalar, float *const out, uint64_t size);
void vec_diff(const float *const vec1, const float *const vec2, float *const out, uint64_t size);
void linear_comb(const float *const vec1, const float *const vec2, const float a1, const float a2, float *const out, uint64_t size);
float l1_norm(const float *const vec, uint64_t size);

void find_dangling(const csr_matrix_t *const mat, uint64_t *const out, uint64_t *out_size);
float rank_loss(const uint64_t *const dangling_indices, const uint64_t dangling_size, const float *const rank, const uint64_t size);
float random_surfer(float *new_rank, float *last_rank, float *e, uint64_t size, float random_jump_prob, float teleport);
int compare_degrees(const void *a, const void *b);
void reorder_by_degree(csr_matrix_t *mat);

int main(int argc, char **argv) {
    // Check errors
    if (argc < 2 || argc > 3) {
        fprintf(stderr, "Usage: %s path/to/web_graph.csr [tries]\n", argv[0]);
        return -1;
    }

    // Parse tries (default to 10 if not provided)
    int tries = 10; 
    if (argc == 3) {
        tries = atoi(argv[2]);
        if (tries <= 0) {
            fprintf(stderr, "Error: Number of tries must be a positive integer.\n");
            return -1;
        }
    }    
    
    // Load web graph
    csr_matrix_t web;
    if (load_csr(argv[1], &web) != 0) {
        printf("failed loading graph\n");
        return -1;
    }
    
#   ifdef GRAPH_REORDERING
    // Graph reordering
    reorder_by_degree(&web);
#   endif

    // Timer data
    execution_times = malloc(sizeof(double) * tries);

    // Preallocate buffers
    uint64_t size = web.n_rows;
    e = malloc(sizeof(float) * size);
    vec0 = malloc(sizeof(float) * size);
    vec1 = malloc(sizeof(float) * size);
    diff = malloc(sizeof(float) * size);
    dangling = malloc(sizeof(uint64_t) * size);

    // Run PageRank
    float *rank = malloc(sizeof(float) * web.n_rows);
    for (int i = 0; i < tries; i++) {
        double start = get_time();
        if (!pagerank(&web, &rank)) {
            cleanup(&web, rank);
            return -1;
        }
        double end = get_time();

        execution_times[i] = end - start;
    }

    // Print metrics
    double avg = get_average(execution_times, tries);
    double med = get_median(execution_times, tries);
    printf("Average time: %lf - Median time: %lf\n", avg, med);

    // Cleanup
    cleanup(&web, rank);

    return 0;
}


/**
 * PageRank implementation
 */
bool pagerank(const csr_matrix_t *const mat, float **rank) {
    if (mat->n_rows != mat->n_columns) {
        fprintf(stderr, "ERROR: n_rows (%" PRIu64 ") must be equal to n_columns (%" PRIu64 ")", mat->n_rows, mat->n_columns);
        return false;
    }

    const float EPSILON = 1e-6;
    const float GLOBAL_RANK = 1.0;
    const uint32_t MAX_ITERATIONS = 1000;
    const float RANDOM_JUMP_PROB = 0.15;
    
    uint64_t size = mat->n_rows;

    // Initial distribution: uniform
    float init_val = GLOBAL_RANK / size;
    for (uint64_t i = 0; i < size; i++) {
        e[i] = init_val;
    }

    // Initialization
    float *last_rank = vec0;
    float *new_rank = vec1;
    memcpy(vec0, e, sizeof(float) * size);
    uint64_t dangling_size;
    find_dangling(mat, dangling, &dangling_size);

    float delta;
    uint32_t iteration = 0;

    do {
        // Matrix-Vector multiplication
        matvec_mul(mat, last_rank, new_rank);
        
        // Compute lost rank (dangling pages)
        float teleport = rank_loss(dangling, dangling_size, new_rank, size);
        
#       ifndef MERGE_VECTOR_OPERATIONS
        // Random surfer
        linear_comb(new_rank, e, 1.0 - RANDOM_JUMP_PROB, RANDOM_JUMP_PROB, new_rank, size);
        vec_add_scalar(new_rank, teleport / (float)size, new_rank, size);

        // Compute delta
        vec_diff(new_rank, last_rank, diff, size);
        delta = l1_norm(diff, size);

#       else
        delta = random_surfer(new_rank, last_rank, e, size, RANDOM_JUMP_PROB, teleport);
#       endif

        // Update data of last iteration
        float *aux = new_rank;
        new_rank = last_rank;
        last_rank = aux;

        iteration++;

        // Check convergence
    } while(delta > EPSILON && iteration < MAX_ITERATIONS);

    // Store result
    memcpy(*rank, last_rank, sizeof(float) * size);

    return true;
}

/**
 * Frees the resources used by the main function
 */
void cleanup(csr_matrix_t *mat, float *rank) {
    free(rank);
    free(mat->values);
    free(mat->columns);
    free(mat->row_ptrs);
    free(e);
    free(vec0);
    free(vec1);
    free(diff);
    free(dangling);
    free(execution_times);
}

int compare_doubles(const void *a, const void *b) {
    double arg1 = *(const double *)a;
    double arg2 = *(const double *)b;
    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

double get_average(double *times, int n) {
    double sum = 0;
    for (int i = 0; i < n; i++) sum += times[i];
    return sum / n;
}

double get_median(double *times, int n) {
    qsort(times, n, sizeof(double), compare_doubles);
    if (n % 2 == 0) {
        return (times[n / 2 - 1] + times[n / 2]) / 2.0;
    } else {
        return times[n / 2];
    }
}



// Helper functions

void matvec_mul(const csr_matrix_t *const mat, const float *const vec, float *const out) {
    #pragma omp parallel for schedule(static, 512)
    for (uint64_t j = 0; j < mat->n_rows; j++) {
        float sum = 0.0;

        for (uint64_t index = mat->row_ptrs[j]; index < mat->row_ptrs[j + 1]; index++) {
            sum += mat->values[index] * vec[mat->columns[index]];
        }

        out[j] = sum;
    }
}

void vec_add_scalar(const float *const vec1, float scalar, float *const out, uint64_t size) {
    #pragma omp parallel for
    for (uint64_t i = 0; i < size; i++) {
        out[i] = vec1[i] + scalar;
    }
}

void vec_diff(const float *const vec1, const float *const vec2, float *const out, uint64_t size) {
    #pragma omp parallel for
    for (uint64_t i = 0; i < size; i++) {
        out[i] = vec1[i] - vec2[i];
    }
}

void linear_comb(const float *const vec1, const float *const vec2, const float a1, const float a2, float *const out, uint64_t size) {
    #pragma omp parallel for
    for (uint64_t i = 0; i < size; i++) {
        out[i] = vec1[i] * a1 + vec2[i] * a2;
    }
}

float l1_norm(const float *const vec, uint64_t size) {
    float sum = 0.0;

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

float rank_loss(const uint64_t *const dangling_indices, const uint64_t dangling_size, const float *const rank, const uint64_t size) {
    float sum = 0.0;

    #pragma omp parallel for reduction(+:sum)
    for (uint64_t i = 0; i < dangling_size; i++) {
        sum += rank[dangling_indices[i]];
    }

    return sum;
}

float random_surfer(float *new_rank, float *last_rank, float *e, uint64_t size, float random_jump_prob, float teleport) {
    float delta = 0.0;

    #pragma omp parallel for reduction(+:delta)
    for (uint64_t i = 0; i < size; i++) {
        // Random Surfer
        float v = (1.0 - random_jump_prob) * new_rank[i] + random_jump_prob * e[i];
        v += teleport / size;
        new_rank[i] = v;

        // Compute Delta
        delta += fabs(v - last_rank[i]);
    }

    return delta;
}

int compare_degrees(const void *a, const void *b) {
    node_degree_t *nodeA = (node_degree_t *)a;
    node_degree_t *nodeB = (node_degree_t *)b;
    
    // Higher degree first
    return (nodeB->degree - nodeA->degree);
}

void reorder_by_degree(csr_matrix_t *mat) {
    uint64_t n = mat->n_rows;
    
    // Resource allocation
    node_degree_t *degrees = malloc(n * sizeof(node_degree_t));
    uint64_t *old_to_new = malloc(n * sizeof(uint64_t));
    
    if (!degrees || !old_to_new) {
        fprintf(stderr, "Allocation failed for reordering metadata\n");
        return;
    }

    // Calculate degrees
    #pragma omp parallel for
    for (uint64_t i = 0; i < n; i++) {
        degrees[i].original_id = i;
        degrees[i].degree = (uint32_t)(mat->row_ptrs[i+1] - mat->row_ptrs[i]);
    }

    // Sort nodes by degree
    qsort(degrees, n, sizeof(node_degree_t), compare_degrees);

    // Create mapping from old index to new index
    #pragma omp parallel for
    for (uint64_t i = 0; i < n; i++) {
        old_to_new[degrees[i].original_id] = i;
    }

    // Allocate new buffers
    uint64_t *new_row_ptrs = malloc((n + 1) * sizeof(uint64_t));
    uint64_t *new_columns = malloc(mat->nnz * sizeof(uint64_t));
    float *new_values = malloc(mat->nnz * sizeof(float));

    if (!new_row_ptrs || !new_columns || !new_values) {
        fprintf(stderr, "Allocation failed for new CSR buffers\n");
        exit(1);
    }

    // Build row_ptrs
    new_row_ptrs[0] = 0;
    for (uint64_t i = 0; i < n; i++) {
        uint64_t old_idx = degrees[i].original_id;
        uint64_t num_edges = mat->row_ptrs[old_idx + 1] - mat->row_ptrs[old_idx];
        new_row_ptrs[i + 1] = new_row_ptrs[i] + num_edges;
    }

    // Fill the data
    
    //#pragma omp parallel for
    for (uint64_t i = 0; i < n; i++) {
        uint64_t old_idx = degrees[i].original_id;
        uint64_t old_start = mat->row_ptrs[old_idx];
        uint64_t new_start = new_row_ptrs[i];
        uint64_t num_edges = mat->row_ptrs[old_idx + 1] - mat->row_ptrs[old_idx];
    
        for (uint64_t j = 0; j < num_edges; j++) {
            uint64_t src_idx = old_start + j;
            uint64_t dst_idx = new_start + j;
            
            // Check source CSR bounds
            if (src_idx >= mat->nnz) {
                printf("Error: src_idx %lu out of bounds (nnz %lu)\n", src_idx, mat->nnz);
                continue; 
            }
    
            uint64_t target_node = mat->columns[src_idx];
    
            // Check mapping bounds (The most likely culprit)
            if (target_node >= n) {
                printf("CRITICAL: target_node %lu >= n %lu. Mapping failed.\n", target_node, n);
                exit(1);
            }
    
            new_columns[dst_idx] = old_to_new[target_node];
            new_values[dst_idx] = mat->values[src_idx];
        }
    }

    // Replace old buffers and cleanup
    free(mat->row_ptrs);
    free(mat->columns);
    free(mat->values);

    mat->row_ptrs = new_row_ptrs;
    mat->columns = new_columns;
    mat->values = new_values;

    free(degrees);
    free(old_to_new);
}
