#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CentHORD: Detect diagonals (HOR) and checkerboard patterns (ancestral dimeric HOR) in distance matrices.


 
Features:
- Diagonal detection with optional merging of nearby fragments.
- Checkerboard detection using normalized cross-correlation (NCC).
- Fusion of overlapping/adjacent checkerboard blocks into larger zones.
- Outputs CSVs for diagonals and checkerboard zones.
- Plots a heatmap with:
    - Non-overlapping diagonals annotated
    - Checkerboard zones highlighted
    - A barplot of NCC scores above the heatmap
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter
from scipy.signal import convolve2d
from scipy.ndimage import uniform_filter
import h5py
import os


# ===============================
# Utility functions
# ===============================

def load_matrix(file_path, input_type):
    """Load a distance matrix from text or HDF5."""
    if input_type == "text":
        return np.loadtxt(file_path).astype(float)
    elif input_type == "h5":
        with h5py.File(file_path, "r") as f:
            dset = list(f.values())[0]  # assumes only one dataset
            return dset[()].astype(float)
    else:
        raise ValueError("Unknown input type. Use 'text' or 'h5'.")
        


def anti_checker_prewhiten(mat, k=2):
    """
    Small k×k mean filter (k must be even) to attenuate checkerboard patterns.
    Keeps scale (normalized kernel). Use k=0 to disable.
    """
    if k <= 0:
        return mat
    ker = np.ones((k, k), dtype=float) / (k * k)
    mat_filled = np.nan_to_num(mat, nan=np.nanmean(mat))
    return convolve2d(mat_filled, ker, mode='same', boundary='symm')


def robust_local_stats(mat, size=7):
    """
    Compute local median and local MAD (converted to SD equivalent).
    Returns (median, sd_equiv).
    """
    mat_filled = np.nan_to_num(mat, nan=np.nanmean(mat))
    med = median_filter(mat_filled, size=size, mode="nearest")
    abs_dev = np.abs(mat_filled - med)
    mad = median_filter(abs_dev, size=size, mode="nearest")
    sd_equiv = 1.4826 * mad
    return med, sd_equiv


def local_stats(mat, size=7):
    """Compute local mean and SD using a square window."""
    mask = ~np.isnan(mat)
    mat_filled = np.nan_to_num(mat, nan=np.nanmean(mat))
    local_mean = uniform_filter(mat_filled, size=size, mode="nearest")
    local_sq_mean = uniform_filter(mat_filled**2, size=size, mode="nearest")
    local_var = np.maximum(local_sq_mean - local_mean**2, 0)
    local_sd = np.sqrt(local_var)
    local_mean[~mask] = np.nan
    local_sd[~mask] = np.nan
    return local_mean, local_sd


def convolve_diagonal(mat, size=5):
    """Apply diagonal convolution with given window size."""
    kernel = np.eye(size)
    mat_filled = np.nan_to_num(mat, nan=np.nanmean(mat))
    return convolve2d(mat_filled, kernel, mode='same', boundary='fill', fillvalue=0)


def merge_segments(diagonals_df, max_gap=5, mode="fixed"):
    """Merge nearby diagonal fragments."""
    if diagonals_df.empty:
        return diagonals_df

    merged = []
    for diag_id, group in diagonals_df.groupby("diag_id"):
        group = group.sort_values(by="start_row").reset_index(drop=True)
        current = group.iloc[0].to_dict()
        for _, row in group.iloc[1:].iterrows():
            gap = int(row["start_row"]) - int(current["end_row"]) - 1
            threshold = max_gap if mode == "fixed" else max_gap * min(current["length"], int(row["length"])) / 100
            if gap <= threshold:
                current["end_row"] = int(row["end_row"])
                current["end_col"] = int(row["end_col"])
                current["length"] = int(row["end_row"]) - int(current["start_row"]) + 1
                current["mid_row"] = int(np.round((current["start_row"] + current["end_row"]) / 2))
                current["mid_col"] = int(np.round((current["start_col"] + current["end_col"]) / 2))
                current["mean_distance"] = np.mean([current["mean_distance"], row["mean_distance"]])
            else:
                merged.append(current)
                current = row.to_dict()
        merged.append(current)
    return pd.DataFrame(merged)


def extract_diagonals(mat, conv_size=5, local_size=7, sd_factor=2.0, min_length=4,
                      merge=False, merge_gap=5, prewhiten_checker=0, local_method="robust", merge_mode="fixed"):
    """Detect significant diagonals based on local mean & SD."""

    # Pre-whitening to attenuate checkerboard patterns
    mat_pw = anti_checker_prewhiten(mat, k=prewhiten_checker)

    # Convolution along diagonals (comme avant, mais sur mat_pw)
    conv_matrix = convolve_diagonal(mat_pw, size=conv_size)

    # Local background estimation
    if local_method == "robust":
        background, bg_sd = robust_local_stats(conv_matrix, size=local_size)
    else:
        background, bg_sd = local_stats(conv_matrix, size=local_size)

    # Relative score
    score = conv_matrix - background
    mask = score < - sd_factor * bg_sd
    coords = np.argwhere(mask)
    diag_ids = coords[:, 1] - coords[:, 0]





    diagonals = []
    for diag_id in np.unique(diag_ids):
        if diag_id <= 0:
            continue
        diag_points = coords[diag_ids == diag_id]
        diag_points = diag_points[np.argsort(diag_points[:, 0])]
        breaks = np.where(np.diff(diag_points[:, 0]) > 1)[0] + 1
        segments = np.split(diag_points, breaks)
        for seg in segments:
            if len(seg) >= min_length:
                vals = [mat[r, c] for r, c in seg]
                diagonals.append({
                    "diag_id": int(diag_id),
                    "start_row": int(seg[0, 0]),
                    "end_row": int(seg[-1, 0]),
                    "start_col": int(seg[0, 1]),
                    "end_col": int(seg[-1, 1]),
                    "length": len(seg),
                    "mid_row": int(np.round(np.mean(seg[:, 0]))),
                    "mid_col": int(np.round(np.mean(seg[:, 1]))),
                    "mean_distance": float(np.mean(vals))
                })
    df = pd.DataFrame(diagonals)
    if merge and not df.empty:
        df = merge_segments(df, max_gap=merge_gap, mode=merge_mode)
    return df


def select_non_overlapping_x_only(diagonals_df):
    """Select the largest non-overlapping diagonals (diag_id > 0 only)."""
    if diagonals_df.empty:
        return pd.DataFrame()
    diagonals_df = diagonals_df[diagonals_df["diag_id"] > 0]
    diagonals_df = diagonals_df.sort_values(by=["length", "diag_id"], ascending=[False, True]).reset_index(drop=True)
    if diagonals_df.empty:
        return pd.DataFrame()
    selected = []
    occupied = np.zeros((int(diagonals_df["end_row"].max()) + 2,), dtype=bool)
    for _, row in diagonals_df.iterrows():
        start = int(row["start_row"])
        end = int(row["end_row"])
        if not np.any(occupied[start:end + 1]):
            selected.append(row)
            occupied[start:end + 1] = True
    return pd.DataFrame(selected)
    
def select_non_overlapping_length_only(diagonals_df):
    """Select the largest non-overlapping diagonals (diag_id > 0 only),
    excluding the main diagonal and negative offsets.
    Non-overlap is enforced on both rows (Y) and columns (X).
    """
    if diagonals_df.empty:
        return pd.DataFrame()
    diagonals_df = diagonals_df[diagonals_df["diag_id"] > 0]
    diagonals_df = diagonals_df.sort_values(by=["length", "diag_id"], ascending=[False, True]).reset_index(drop=True)
    if diagonals_df.empty:
        return pd.DataFrame()

    selected = []
    # Track occupancy in both row and col axes
    max_row = int(diagonals_df["end_row"].max()) + 2
    max_col = int(diagonals_df["end_col"].max()) + 2
    occupied_rows = np.zeros(max_row, dtype=bool)
    occupied_cols = np.zeros(max_col, dtype=bool)

    for _, row in diagonals_df.iterrows():
        r_start, r_end = int(row["start_row"]), int(row["end_row"])
        c_start, c_end = int(row["start_col"]), int(row["end_col"])

        if (not np.any(occupied_rows[r_start:r_end + 1]) and
            not np.any(occupied_cols[c_start:c_end + 1])):
            selected.append(row)
            occupied_rows[r_start:r_end + 1] = True
            occupied_cols[c_start:c_end + 1] = True

    return pd.DataFrame(selected)
    
def select_non_overlapping(diagonals_df, criterion="length"):
    """Select non-overlapping diagonals based on criterion.
    criterion: "length" (default) or "distance" (mean distance).
    """
    if diagonals_df.empty:
        return pd.DataFrame()
    diagonals_df = diagonals_df[diagonals_df["diag_id"] > 0]

    if criterion == "length":
        diagonals_df = diagonals_df.sort_values(
            by=["length", "diag_id"], ascending=[False, True]
        ).reset_index(drop=True)
    elif criterion == "distance":
        diagonals_df = diagonals_df.sort_values(
            by=["mean_distance", "diag_id"], ascending=[True, True]
        ).reset_index(drop=True)

    if diagonals_df.empty:
        return pd.DataFrame()

    selected = []
    max_row = int(diagonals_df["end_row"].max()) + 2
    max_col = int(diagonals_df["end_col"].max()) + 2
    occupied_rows = np.zeros(max_row, dtype=bool)
    occupied_cols = np.zeros(max_col, dtype=bool)

    for _, row in diagonals_df.iterrows():
        r_start, r_end = int(row["start_row"]), int(row["end_row"])
        c_start, c_end = int(row["start_col"]), int(row["end_col"])

        if (not np.any(occupied_rows[r_start:r_end + 1]) and
            not np.any(occupied_cols[c_start:c_end + 1])):
            selected.append(row)
            occupied_rows[r_start:r_end + 1] = True
            occupied_cols[c_start:c_end + 1] = True

    return pd.DataFrame(selected)


# ===============================
# Checkerboard detection (NCC)
# ===============================

def make_checker_kernelOLD(size=6):
    """Generate checkerboard kernel (upper-triangular, zero mean)."""
    if size % 2 != 0:
        raise ValueError("Checkerboard kernel size must be even.")
    block = np.array([[1, -1], [-1, 1]])
    reps = size // 2
    kernel = np.tile(block, (reps, reps))
    kernel = np.triu(kernel, 1)  # keep only upper triangle (excluding diagonal)
    kernel = kernel - np.mean(kernel[kernel != 0])
    return kernel

def make_checker_kernel(size=7):
    """Generate a checkerboard kernel of given odd size."""
    if size % 2 == 0:
        raise ValueError("Checkerboard kernel size must be odd.")
    pattern = np.indices((size, size)).sum(axis=0) % 2
    kernel = pattern * 2 - 1  # values in {-1, 1}
    return kernel

def detect_checkerboardOLD(mat, checker_size=6, ncc_threshold=0.5):
    """Detect checkerboard patterns with NCC (values in [-1, 1])."""
    # Ensure kernel size is odd
    if checker_size % 2 == 0:
        checker_size += 1

    kernel = make_checker_kernel(checker_size)
    mat_filled = np.nan_to_num(mat, nan=np.nanmean(mat))

    nrows, ncols = mat_filled.shape
    ncc = np.full_like(mat_filled, np.nan, dtype=float)
    kh, kw = kernel.shape
    kh2, kw2 = kh // 2, kw // 2
    norm_kernel = np.linalg.norm(kernel)

    for i in range(kh2, nrows - kh2):
        for j in range(kw2, ncols - kw2):
            patch = mat_filled[i - kh2:i + kh2 + 1, j - kw2:j + kw2 + 1]
            if patch.shape != kernel.shape:
                continue
            norm_patch = np.linalg.norm(patch)
            if norm_patch > 0:
                ncc[i, j] = np.sum(patch * kernel) / (norm_patch * norm_kernel)

    hits = np.argwhere(ncc > ncc_threshold)
    return ncc, hits

def detect_checkerboardOLD(mat, checker_size=7, ncc_threshold=0.5, band_size=0):
    """
    Detect checkerboard patterns with NCC along the main diagonal only.
    - checker_size: size of the checkerboard kernel (must be odd)
    - ncc_threshold: threshold for NCC score
    - band_size: number of rows/cols around the main diagonal to include (0 = exactly diagonal)
    """
     # Ensure kernel size is odd
    if checker_size % 2 == 0:
        checker_size += 1

    kernel = make_checker_kernel(checker_size)
    mat_filled = np.nan_to_num(mat, nan=np.nanmean(mat))

    nrows, ncols = mat_filled.shape
    ncc = np.full_like(mat_filled, np.nan, dtype=float)
    kh, kw = kernel.shape
    kh2, kw2 = kh // 2, kw // 2
    norm_kernel = np.linalg.norm(kernel)

    # Scan only along the main diagonal (and optional band around it)
    for i in range(kh2, nrows - kh2):
        j = i  # position along the diagonal
        for offset in range(-band_size, band_size + 1):
            jj = j + offset
            if jj - kw2 < 0 or jj + kw2 + 1 > ncols:
                continue
            patch = mat_filled[i - kh2:i + kh2 + 1, jj - kw2:jj + kw2 + 1]
            if patch.shape != kernel.shape:
                continue
            norm_patch = np.linalg.norm(patch)
            if norm_patch > 0:
                ncc[i, jj] = np.sum(patch * kernel) / (norm_patch * norm_kernel)

    hits = np.argwhere(ncc > ncc_threshold)
    return ncc, hits

def merge_checker_hitsOLD(hits, max_gap=5):
    """Merge overlapping/adjacent checkerboard hits."""
    if len(hits) == 0:
        return []
    hits = sorted(hits, key=lambda x: (x[0], x[1]))
    merged = []
    current = [hits[0][0], hits[0][1], hits[0][0], hits[0][1]]
    for (i, j) in hits[1:]:
        if i <= current[2] + max_gap and j <= current[3] + max_gap:
            current[2] = max(current[2], i)
            current[3] = max(current[3], j)
        else:
            merged.append(tuple(current))
            current = [i, j, i, j]
    merged.append(tuple(current))
    return merged



def detect_checkerboard(mat, checker_size=7, ncc_threshold=0.5, band_size=0):
    """
    Detect checkerboard patterns with NCC along the main diagonal only.
    Returns both the NCC matrix and detected zones (with kernel size included).
    """
    if checker_size % 2 == 0:
        checker_size += 1

    kernel = make_checker_kernel(checker_size)
    mat_filled = np.nan_to_num(mat, nan=np.nanmean(mat))

    nrows, ncols = mat_filled.shape
    ncc = np.full_like(mat_filled, np.nan, dtype=float)
    kh, kw = kernel.shape
    kh2, kw2 = kh // 2, kw // 2
    norm_kernel = np.linalg.norm(kernel)

    zones = []
    for i in range(kh2, nrows - kh2):
        j = i
        for offset in range(-band_size, band_size + 1):
            jj = j + offset
            if jj - kw2 < 0 or jj + kw2 + 1 > ncols:
                continue
            patch = mat_filled[i - kh2:i + kh2 + 1, jj - kw2:jj + kw2 + 1]
            if patch.shape != kernel.shape:
                continue
            norm_patch = np.linalg.norm(patch)
            if norm_patch > 0:
                ncc_val = np.sum(patch * kernel) / (norm_patch * norm_kernel)
                ncc[i, jj] = ncc_val
                if ncc_val > ncc_threshold:
                    # Include kernel borders
                    r1, r2 = i - kh2, i + kh2
                    c1, c2 = jj - kw2, jj + kw2
                    zones.append((r1, c1, r2, c2))

    return ncc, zones


def merge_checker_hits(zones, max_gap=5):
    """Merge overlapping/adjacent checkerboard zones (rectangles)."""
    if len(zones) == 0:
        return []

    zones = sorted(zones, key=lambda x: (x[0], x[1]))
    merged = []
    current = list(zones[0])
    for (r1, c1, r2, c2) in zones[1:]:
        if r1 <= current[2] + max_gap and c1 <= current[3] + max_gap:
            current[2] = max(current[2], r2)
            current[3] = max(current[3], c2)
        else:
            merged.append(tuple(current))
            current = [r1, c1, r2, c2]
    merged.append(tuple(current))
    return merged





# ===============================
# Plotting
# ===============================

def plot_heatmap(dist_matrix, nonoverlapping_df, diagonals_df, checker_zones, ncc_scores,
                 output_png, invert_colors=False, min_diag_len=4, diag_offset=1, diag_width=3, diag_color="red"):
    """Plot heatmap + diagonals + checkerboard zones + NCC barplot."""
    cmap_name = "RdYlBu_r" if not invert_colors else "YlOrRd"
        
    cmap = plt.colormaps[cmap_name]
    fig = plt.figure(figsize=(10, 12), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 10], hspace=0.05)

    # NCC barplot
    ax_bar = fig.add_subplot(gs[0])
    if ncc_scores is not None:
        x = np.arange(len(ncc_scores))
        ax_bar.bar(x, ncc_scores, color="blue", width=1.0)
        ax_bar.set_xlim(0, dist_matrix.shape[1])
        ax_bar.set_ylabel("NCC")
        ax_bar.set_xticks([])
    else:
        ax_bar.axis("off")

    # Heatmap
    ax = fig.add_subplot(gs[1])
    im = ax.imshow(dist_matrix, origin="upper", cmap=cmap)
    fig.colorbar(im, ax=ax, label="Distance")
    # Diagonals
    
    if not nonoverlapping_df.empty:
        for _, row in nonoverlapping_df.iterrows():
            if row["length"] >= min_diag_len:
                ax.plot([row["start_col"] + diag_offset, row["end_col"] + diag_offset],
                        [row["start_row"], row["end_row"]],
                        color=diag_color, lw=diag_width)
                ax.text(row["mid_col"] + diag_offset, row["mid_row"], str(int(row["length"])),
                        color="black", fontsize=8, ha="center", va="center",
                        bbox=dict(facecolor="white", alpha=0.6, edgecolor="none", pad=1))

    for _, row in diagonals_df.iterrows():
        x = np.arange(row["start_col"], row["end_col"] + 1) + diag_offset
        y = np.arange(row["start_row"], row["end_row"] + 1)
        # keep only points in lower-left triangle (y >= x)
        #print(x)
        mask = x >= y
        ax.plot(y[mask], x[mask],
             color=diag_color,
             linewidth=diag_width)
             
             
    # Checkerboard zones
    for (r1, c1, r2, c2) in checker_zones:
        rect = plt.Rectangle((c1, r1), c2 - c1, r2 - r1,
                             edgecolor="blue", facecolor="none", lw=2)
        ax.add_patch(rect)

    ax.set_xlabel("Columns")
    ax.set_ylabel("Rows")
    ax.set_title("Distance matrix with diagonals and checkerboard zones")
 #   plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    plt.close()


# ===============================
# Main
# ===============================

def main():
    parser = argparse.ArgumentParser(description="Detect diagonals and checkerboards in a distance matrix")
    parser.add_argument("--input", required=True, help="Input distance matrix file (.txt or .h5)")
    parser.add_argument("--input_type", choices=["text", "h5"], required=True, help="Input file type")
    parser.add_argument("--invert_colors", action="store_true", help="Invert heatmap colors")
    parser.add_argument("--detect_diagonals", action="store_true", help="Enable diagonal detection")
    parser.add_argument("--conv_size", type=int, default=5, help="Diagonal convolution kernel size (odd)")
    parser.add_argument("--local_size", type=int, default=7, help="Local window size for diagonal stats")
    parser.add_argument("--sd_factor", type=float, default=2.0, help="SD factor for diagonal detection")
    parser.add_argument("--min_diag_len", type=int, default=4, help="Minimum diagonal length to annotate")
    parser.add_argument("--merge_diagonals", action="store_true", help="Enable merging of diagonal fragments")
    parser.add_argument("--merge_gap", type=int, default=5, help="Gap for merging diagonals")
    parser.add_argument("--merge_mode", choices=["fixed", "relative"], default="fixed", help="Merge mode")
    parser.add_argument("--local_method", choices=["classic", "robust"], default="robust",
                    help="Method for local background estimation (default: robust = median+MAD).")
    parser.add_argument("--prewhiten_checker", type=int, default=0,
                    help="Apply prewhitening (small mean filter) to reduce checkerboard effect. "
                         "Set to 2 or 4 for strong checkerboard. Default=0 (off).")
    parser.add_argument("--diag_select", choices=["length", "distance"], default="length",
                    help="Criterion for selecting non-overlapping diagonals: "
                         "'length' = longest diagonals first (default), "
                         "'distance' = smallest mean distance first.") 
    parser.add_argument("--diag_offset", type=int, default=0,
                    help="Horizontal offset applied when plotting diagonals")
    parser.add_argument("--diag_color", type=str, default="black",
                    help="Color of diagonal lines on the heatmap")
    parser.add_argument("--diag_width", type=float, default=1.5,
                    help="Line width of diagonal lines on the heatmap")

    # Checkerboard options
    parser.add_argument("--detect_checkerboard", action="store_true", help="Enable checkerboard detection")
    parser.add_argument("--checker_size", type=int, default=11, help="Checkerboard kernel size (odd)")
    parser.add_argument("--checker_ncc_threshold", type=float, default=0.5, help="NCC threshold for checkerboards")
    parser.add_argument("--checker_merge_gap", type=int, default=5, help="Gap for merging checkerboard zones")

    args = parser.parse_args()
    if args.checker_size % 2 == 0:
        print("Error: checkerboard kernel size must be odd (size=%i)" % (args.checker_size))
        quit()
        
    dist_matrix = load_matrix(args.input, args.input_type)
    prefix = os.path.splitext(os.path.basename(args.input))[0]

    diagonals_df, nonoverlapping_df = pd.DataFrame(), pd.DataFrame()
    checker_zones, ncc_scores = [], None

    if args.detect_diagonals:
        diagonals_df = extract_diagonals(dist_matrix,
                                         conv_size=args.conv_size,
                                         local_size=args.local_size,
                                         sd_factor=args.sd_factor,
                                         min_length=args.min_diag_len,
                                         merge=args.merge_diagonals,
                                         merge_gap=args.merge_gap,
                                         prewhiten_checker=args.prewhiten_checker, 
                                         local_method=args.local_method,
                                         merge_mode=args.merge_mode)
        diagonals_df.to_csv(prefix + "_diagonals.csv", index=False)
        nonoverlapping_df = select_non_overlapping(diagonals_df, args.diag_select)
        nonoverlapping_df.to_csv(prefix + "_nonoverlapping_diagonals.csv", index=False)
        
        
    if args.detect_checkerboard:
            ncc, hits = detect_checkerboard(dist_matrix,
                                    checker_size=args.checker_size,
                                    ncc_threshold=args.checker_ncc_threshold)
            checker_zones = merge_checker_hits(hits, max_gap=args.checker_merge_gap)

            # Extract NCC scores along main diagonal
            ncc_diag = np.array([ncc[i, i] for i in range(min(ncc.shape))])
            ncc_scores = np.pad(ncc_diag, (0, dist_matrix.shape[0] - len(ncc_diag)), constant_values=np.nan)

            # Compute mean NCC per zone  ### NEW
            results = []
            for (r1, c1, r2, c2) in checker_zones:
                sub_ncc = ncc[r1:r2+1, c1:c2+1]
                mean_ncc = float(np.nanmean(sub_ncc))
                results.append([r1, c1, r2, c2, mean_ncc])

            pd.DataFrame(results, columns=["r1", "c1", "r2", "c2", "mean_ncc"]).to_csv(
                prefix + "_checkerboard_zones.csv", index=False)
        
        
    plot_heatmap(dist_matrix, nonoverlapping_df, diagonals_df, checker_zones, ncc_scores,
                  prefix + "_heatmap.png", invert_colors=args.invert_colors,
                  min_diag_len=args.min_diag_len, diag_offset=args.diag_offset, diag_color=args.diag_color, diag_width=args.diag_width)

    # if args.detect_checkerboard:
        # ncc, hits = detect_checkerboard(dist_matrix,
                                        # checker_size=args.checker_size,
                                        # ncc_threshold=args.checker_ncc_threshold)
        # checker_zones = merge_checker_hits(hits, max_gap=args.checker_merge_gap)
        # ncc_diag = np.array([ncc[i, i] for i in range(min(ncc.shape))])
        # ncc_scores = np.pad(ncc_diag, (0, dist_matrix.shape[0] - len(ncc_diag)), constant_values=np.nan)
        # pd.DataFrame(checker_zones, columns=["r1", "c1", "r2", "c2"]).to_csv(
            # prefix + "_checkerboard_zones.csv", index=False)



if __name__ == "__main__":
    main()
