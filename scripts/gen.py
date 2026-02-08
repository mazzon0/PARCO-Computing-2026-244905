import numpy as np
import struct
import argparse
import os
from scipy.sparse import csr_matrix

def create_graph(n_nodes, n_nnz):
    # Generate edges as (row, col) pairs directly
    # Using a set to ensure unique edges if NNZ is large relative to N
    # For very sparse graphs (web graphs), duplicates are rare enough 
    # that we can just generate and unique them.
    
    print(f"Sampling {n_nnz} edges...")
    
    # Generate random endpoints
    rows = np.random.randint(0, n_nodes, size=n_nnz, dtype=np.uint64)
    cols = np.random.randint(0, n_nodes, size=n_nnz, dtype=np.uint64)
    
    # PageRank doesn't usually care about self-loops, but if you want to
    # remove them or handle duplicates, we do it here:
    edges = np.unique(np.stack((rows, cols), axis=1), axis=0)
    
    actual_nnz = edges.shape[0]
    if actual_nnz < n_nnz:
        print(f"Note: Removed {n_nnz - actual_nnz} duplicate edges.")
    
    data = np.ones(actual_nnz, dtype=np.float32)
    
    print("Converting to CSR format...")
    return csr_matrix((data, (edges[:, 0], edges[:, 1])), shape=(n_nodes, n_nodes))

def save_as_binary_csr(matrix, base_filename):
    filename = f"{base_filename}.csr"
    n_rows, n_cols = matrix.shape
    nnz = matrix.nnz
    
    # Accessing underlying numpy arrays of the CSR matrix
    with open(filename, 'wb') as f:
        f.write(struct.pack('<QQQ', n_rows, n_cols, nnz))
        f.write(matrix.data.astype(np.float32).tobytes())
        f.write(matrix.indices.astype(np.uint64).tobytes())
        f.write(matrix.indptr.astype(np.uint64).tobytes())
    print(f"Successfully wrote {filename}")

def save_as_edge_list(matrix, base_filename):
    filename = f"{base_filename}.el"
    rows, cols = matrix.nonzero()
    # Write in chunks to save memory during file I/O
    with open(filename, 'w') as f:
        for i in range(len(rows)):
            f.write(f"{rows[i]} {cols[i]}\n")
    print(f"Successfully wrote {filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--nodes", type=int, required=True)
    parser.add_argument("-z", "--nnz", type=int, required=True)
    parser.add_argument("-f", "--filename", type=str, required=True)
    args = parser.parse_args()

    base_name = os.path.splitext(args.filename)[0]
    sparse_matrix = create_graph(args.nodes, args.nnz)
    
    save_as_binary_csr(sparse_matrix, base_name)
    save_as_edge_list(sparse_matrix, base_name)