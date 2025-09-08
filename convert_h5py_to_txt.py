#!/usr/bin/env python3

import argparse
import h5py
import numpy as np

def h5_to_text(h5_file, out_file, dataset_name):
    """
    Lit un dataset dans un fichier HDF5 et l'écrit en matrice carrée texte.
    """
    with h5py.File(h5_file, "r") as h5f:
        if dataset_name not in h5f:
            raise ValueError(f"Dataset '{dataset_name}' not found in {h5_file}")
        mat = h5f[dataset_name][:]
    
    # Écriture texte tabulée
    np.savetxt(out_file, mat, fmt="%d", delimiter="\t")
    print(f"Matrix {mat.shape[0]}x{mat.shape[1]} written in {out_file}")

def main():
    parser = argparse.ArgumentParser(description="Convertir une matrice HDF5 en texte")
    parser.add_argument("-i", "--input", required=True, help="Fichier HDF5 d'entrée")
    parser.add_argument("-o", "--output", required=True, help="Fichier texte de sortie")
    parser.add_argument("-d", "--dataset", default="dist_matrix",
                        help="Name of the file in HDF5 format (default: dist_matrix)")
    
    args = parser.parse_args()
    h5_to_text(args.input, args.output, args.dataset)

if __name__ == "__main__":
    main()
