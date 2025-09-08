#!/usr/bin/env python3
"""
compute_block_matrices.py
Compute dense distance matrices for sequences grouped by contig.
If the number of sequences in a contig exceeds a threshold, split them into
overlapping blocks of size L (block_length) with a step size S (step_size).
Each block is saved as a full dense matrix in its own HDF5 file.
"""

import argparse
import h5py
import numpy as np
from Bio import SeqIO
from collections import defaultdict
import Levenshtein

def read_fasta_by_contig(fasta_path, sep=" "):
    """
    Read a FASTA file and group sequences by contig.
    Each FASTA header is expected to have:
        >id contig_name position
    Returns: {contig: [(position, sequence), ...]}
    """
    contigs = defaultdict(list)
    for record in SeqIO.parse(fasta_path, "fasta"):
        header_parts = record.description.split(sep)
        if len(header_parts) < 3:
            raise ValueError(f"Invalid header: {record.description}")
        contig = header_parts[1]
        position = int(header_parts[2])
        contigs[contig].append((position, str(record.seq)))
    return contigs

def compute_distance_matrix(sequences):
    """
    Compute a dense symmetric distance matrix for the given list of sequences.
    """
    n = len(sequences)
    matrix = np.zeros((n, n), dtype=np.uint16)
    for i in range(n):
        for j in range(i+1, n):
            dist = Levenshtein.distance(sequences[i], sequences[j])
            matrix[i, j] = dist
            matrix[j, i] = dist
    return matrix

def save_matrix_h5(matrix, out_file):
    """
    Save a dense distance matrix to an HDF5 file.
    """
    with h5py.File(out_file, "w") as h5f:
        h5f.create_dataset("dist_matrix", data=matrix, compression="gzip")

def process_contig_full(contig_name, entries, out_prefix, threshold, block_length, step_size):
    """
    Process a contig:
    - If num_sequences <= threshold: compute one full matrix for the whole contig.
    - If num_sequences > threshold: split into blocks of length L with step S and compute each block separately.
    """
    # Sort sequences by position in contig
    entries.sort(key=lambda x: x[0])
    seqs = [seq for _, seq in entries]
    n = len(seqs)

    if n <= threshold:
        print(f"{contig_name}: {n} sequences → full matrix")
        matrix = compute_distance_matrix(seqs)
        out_file = f"{out_prefix}_{contig_name}.h5"
        save_matrix_h5(matrix, out_file)
    else:
        print(f"{contig_name}: {n} sequences → split into blocks of {block_length} with step {step_size}")
        block_idx = 0
        for start in range(0, n - block_length + 1, step_size):
            block_seqs = seqs[start:start + block_length]
            matrix = compute_distance_matrix(block_seqs)
            out_file = f"{out_prefix}_{contig_name}_block{block_idx}.h5"
            save_matrix_h5(matrix, out_file)
            print(f"  Block {block_idx}: sequences {start}–{start + block_length - 1} saved to {out_file}")
            block_idx += 1

def main():
    parser = argparse.ArgumentParser(description="Compute full distance matrices per contig, with block splitting for large contigs.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--out_prefix", required=True, help="Output prefix for HDF5 files")
    parser.add_argument("-t", "--threshold", type=int, default=500, help="Max sequences before splitting into blocks (default: 500)")
    parser.add_argument("-l", "--block_length", type=int, default=500, help="Number of sequences per block (default: 500)")
    parser.add_argument("-s", "--step_size", type=int, default=250, help="Step size between blocks (default: 250)")
    parser.add_argument("--sep", default=" ", help="Separator in FASTA headers (default: space)")
    args = parser.parse_args()

    contigs = read_fasta_by_contig(args.input, sep=args.sep)
    for contig, entries in contigs.items():
        process_contig_full(contig, entries, args.out_prefix, args.threshold, args.block_length, args.step_size)

if __name__ == "__main__":
    main()
