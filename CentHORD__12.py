#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CentHORD: Detect diagonals (HOR) in distance matrices + checkerboard by parity-trim.

Outputs per input:
  <prefix>_diagonals.tsv
  <prefix>_nonoverlapping_diagonals.tsv
  <prefix>_checkerboard_zones.tsv
  <prefix>_heatmap.png

Optional aggregated (when --combined_outputs is set):
  combined_diagonals.tsv
  combined_nonoverlapping_diagonals.tsv
  combined_checkerboard_zones.tsv
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter, generic_filter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
import h5py
import os
from concurrent.futures import ProcessPoolExecutor, as_completed


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


def robust_local_stats(mat, size=7):
    """Compute local median and local MAD (converted to SD equivalent)."""
    mat_filled = np.nan_to_num(mat, nan=np.nanmean(mat))
    med = median_filter(mat_filled, size=size, mode="nearest")
    abs_dev = np.abs(mat_filled - med)
    mad = median_filter(abs_dev, size=size, mode="nearest")
    sd_equiv = 1.4826 * mad
    return med, sd_equiv


def convolve_diagonal(mat, size=5):
    """Apply robust diagonal convolution using the median of diagonal elements."""
    fp = np.eye(size, dtype=bool)
    mat_filled = np.nan_to_num(mat, nan=np.nanmean(mat))
    return generic_filter(mat_filled, function=np.nanmedian, footprint=fp, mode="nearest")


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
                      merge=False, merge_gap=5, merge_mode="fixed"):
    """Detect significant diagonals based on local robust statistics."""
    conv_matrix = convolve_diagonal(mat, size=conv_size)
    background, bg_sd = robust_local_stats(conv_matrix, size=local_size)

    score = conv_matrix - background
    mask = score < - sd_factor * bg_sd
    coords = np.argwhere(mask)
    if coords.size == 0:
        return pd.DataFrame(columns=[
            "diag_id","start_row","end_row","start_col","end_col","length","mid_row","mid_col","mean_distance"
        ])
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


def select_non_overlapping(diagonals_df, criterion="repet"):
    """
    Select non-overlapping diagonals based on a criterion:
      - 'repet'   : maximize length / diag_id (default)
      - 'length'  : longest first
      - 'distance': smallest mean_distance first
    Non-overlap is enforced on both rows and columns.
    """
    if diagonals_df.empty:
        return pd.DataFrame()

    diagonals_df = diagonals_df[diagonals_df["diag_id"] > 0].copy()
    if diagonals_df.empty:
        return pd.DataFrame()

    if criterion == "repet":
        tmp = diagonals_df.assign(score=diagonals_df["length"] / diagonals_df["diag_id"])
        diagonals_df = tmp.sort_values(
            by=["score", "length", "diag_id"], ascending=[False, False, True]
        ).drop(columns=["score"]).reset_index(drop=True)
    elif criterion == "length":
        diagonals_df = diagonals_df.sort_values(
            by=["length", "diag_id"], ascending=[False, True]
        ).reset_index(drop=True)
    elif criterion == "distance":
        diagonals_df = diagonals_df.sort_values(
            by=["mean_distance", "diag_id"], ascending=[True, True]
        ).reset_index(drop=True)

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


def compute_region_summaries(dist_matrix, row):
    """Compute mean and median values for full, off-diagonal, and diagonal regions."""
    r1, r2 = int(row["start_row"]), int(row["end_row"])
    c1, c2 = int(row["start_col"]), int(row["end_col"])

    nrows, ncols = dist_matrix.shape
    r1 = max(0, min(r1, nrows - 1))
    r2 = max(0, min(r2, nrows - 1))
    c1 = max(0, min(c1, ncols - 1))
    c2 = max(0, min(c2, ncols - 1))

    if r2 < r1 or c2 < c1:
        return (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

    region = dist_matrix[r1:r2+1, c1:c2+1]

    vals_full = region.ravel()
    vals_full = vals_full[~np.isnan(vals_full)]
    region_mean_full = np.nan if vals_full.size == 0 else float(np.mean(vals_full))
    region_median_full = np.nan if vals_full.size == 0 else float(np.median(vals_full))

    mask = np.ones(region.shape, dtype=bool)
    np.fill_diagonal(mask, False)
    vals_excl = region[mask]
    vals_excl = vals_excl[~np.isnan(vals_excl)]
    region_mean_excl = np.nan if vals_excl.size == 0 else float(np.mean(vals_excl))
    region_median_excl = np.nan if vals_excl.size == 0 else float(np.median(vals_excl))

    diag_vals = np.diag(region)
    diag_vals = diag_vals[~np.isnan(diag_vals)]
    diagline_mean = np.nan if diag_vals.size == 0 else float(np.mean(diag_vals))
    diagline_median = np.nan if diag_vals.size == 0 else float(np.median(diag_vals))

    return (region_mean_excl, region_median_excl,
            region_mean_full, region_median_full,
            diagline_mean, diagline_median)


# ===============================
# Peeling helpers
# ===============================

def _bg_stat(vals: np.ndarray, combine="median"):
    """Combine background: 'median'|'mean'|'pXX'|'qXX'."""
    if vals.size == 0 or np.all(np.isnan(vals)):
        return np.nan
    c = str(combine).lower()
    if c == "median":
        return np.nanmedian(vals)
    if c == "mean":
        return np.nanmean(vals)
    if c.startswith("p"):
        pct = float(c[1:])
        pct = np.clip(pct, 0, 100)
        return np.nanpercentile(vals, pct)
    if c.startswith("q"):
        q = float(c[1:])
        q = np.clip(q, 0.0, 1.0)
        return np.nanquantile(vals, q)
    raise ValueError(f"combine='{combine}' not recognized.")


def extract_main_diag_segments(mask: np.ndarray, side="upper", min_len=1):
    """
    Extract contiguous segments along j-i=const within the chosen triangle.
    Returns list of Nx2 integer arrays of (row, col).
    """
    n, m = mask.shape
    assert n == m, "Mask must be square."
    if side == "upper":
        tri = np.triu(np.ones((n, n), dtype=bool), k=0)
    else:
        tri = np.tril(np.ones((n, n), dtype=bool), k=0)
    M = mask.astype(bool) & tri
    coords = np.argwhere(M)
    if coords.size == 0:
        return []
    deltas = coords[:, 1] - coords[:, 0]
    segs = []
    for d in np.unique(deltas):
        if side == "upper" and d < 0:
            continue
        if side == "lower" and d > 0:
            continue
        pts = coords[deltas == d]
        pts = pts[np.argsort(pts[:, 0])]
        breaks = np.where(np.diff(pts[:, 0]) > 1)[0] + 1
        for chunk in np.split(pts, breaks):
            if len(chunk) >= min_len:
                segs.append(chunk)
    return segs


def peel_mask_by_contrast(
    dist: np.ndarray,
    mask: np.ndarray,
    side="upper",
    k_offsets=(1, 2),
    min_len=5,
    max_peel_per_side=None,
    combine="p10",
):
    """
    Iterative end-peeling: remove ends while dist[end] > BG(parallel bands ±k).
    """
    n = dist.shape[0]
    U = np.triu(np.ones((n, n), dtype=bool), k=0)
    L = np.tril(np.ones((n, n), dtype=bool), k=0)
    tri = U if side == "upper" else L

    m = (mask > 0).astype(np.uint8).copy()
    m[~tri] = 0

    if isinstance(k_offsets, int):
        k_list = [k_offsets]
    else:
        k_list = list(k_offsets)

    def bg_vals_for_delta(delta: int, i_min: int, i_max: int):
        vals = []
        for k in k_list:
            for d_alt in (delta - k, delta + k):
                ii0 = max(i_min, max(0, -d_alt))
                ii1 = min(i_max, min(n - 1, n - 1 - d_alt))
                if ii1 >= ii0:
                    jj0 = ii0 + d_alt
                    jj1 = ii1 + d_alt
                    vals.append(dist[ii0:ii1 + 1, jj0:jj1 + 1].diagonal(0))
        if not vals:
            return np.array([], dtype=float)
        return np.concatenate(vals)

    changed = True
    while changed:
        changed = False
        segs = extract_main_diag_segments(m, side=side, min_len=min_len)
        if not segs:
            break
        new_m = m.copy()
        for seg in segs:
            ii = seg[:, 0].astype(int)
            jj = seg[:, 1].astype(int)

            peeled_head = 0
            peeled_tail = 0

            while len(ii) >= min_len:
                delta = int(jj[0] - ii[0])
                i_min, i_max = int(ii.min()), int(ii.max())
                vals_bg = bg_vals_for_delta(delta, i_min, i_max)
                stat_bg = _bg_stat(vals_bg, combine=combine)
                if not np.isfinite(stat_bg):
                    break

                did_peel = False
                if (dist[ii[0], jj[0]] > stat_bg) and (max_peel_per_side is None or peeled_head < max_peel_per_side):
                    new_m[ii[0], jj[0]] = 0
                    ii = ii[1:]; jj = jj[1:]
                    peeled_head += 1
                    did_peel = True

                if len(ii) >= min_len:
                    if (dist[ii[-1], jj[-1]] > stat_bg) and (max_peel_per_side is None or peeled_tail < max_peel_per_side):
                        new_m[ii[-1], jj[-1]] = 0
                        ii = ii[:-1]; jj = jj[:-1]
                        peeled_tail += 1
                        did_peel = True

                if not did_peel:
                    break

        if not np.array_equal(new_m, m):
            m = new_m
            changed = True

    # purge final short segments
    for seg in extract_main_diag_segments(m, side=side, min_len=1):
        if len(seg) < min_len:
            m[seg[:, 0], seg[:, 1]] = 0

    segs_after = extract_main_diag_segments(m, side=side, min_len=min_len)
    return m, segs_after


# ===============================
# Checkerboard (parity + trim)
# ===============================

def band_vals_in_square(mat, k, mid, half):
    """
    Values along band offset +k inside square centered at (mid,mid) with half-size 'half'.
    """
    H, W = mat.shape
    r1 = max(0, mid - half)
    r2 = min(H - 1, mid + half)
    c1 = r1
    c2 = r2
    r_start = max(r1, r1 - k)
    r_end   = min(r2, r2 - k)
    if r_end < r_start:
        return np.array([], dtype=float)
    rows = np.arange(r_start, r_end + 1, dtype=int)
    cols = rows + k
    vals = mat[rows, cols]
    return vals.astype(float)


def trim_top(vals, frac):
    if vals.size == 0 or frac <= 0:
        return vals
    q = np.nanquantile(vals, 1.0 - frac)
    return vals[vals <= q]


def trim_bottom(vals, frac):
    if vals.size == 0 or frac <= 0:
        return vals
    q = np.nanquantile(vals, frac)
    return vals[vals >= q]


def merge_overlapping_or_touching(boxes):
    """Merge rectangles that overlap OR touch (adjacent) in either axis."""
    if not boxes:
        return []
    boxes = [list(b) for b in boxes]
    used = [False]*len(boxes)
    merged = []

    def overlap_or_touch(b1, b2):
        return not (b1[2] < b2[0] - 1 or b1[0] > b2[2] + 1 or
                    b1[3] < b2[1] - 1 or b1[1] > b2[3] + 1)

    for i in range(len(boxes)):
        if used[i]:
            continue
        r1, c1, r2, c2 = boxes[i]
        changed = True
        while changed:
            changed = False
            for j in range(i+1, len(boxes)):
                if used[j]:
                    continue
                if overlap_or_touch((r1,c1,r2,c2), boxes[j]):
                    r1 = min(r1, boxes[j][0])
                    c1 = min(c1, boxes[j][1])
                    r2 = max(r2, boxes[j][2])
                    c2 = max(c2, boxes[j][3])
                    used[j] = True
                    changed = True
        used[i] = True
        merged.append((r1, c1, r2, c2))
    return merged


def detect_checkerboard_parity_trimmed_sliding(
    mat,
    square_size=11,       # forced odd; bands are fixed: even_k=2, odd_k=1
    trim_even_top=0.20,   # drop top 20% of even band
    trim_odd_bottom=0.20  # drop bottom 20% of odd band
):
    X = np.nan_to_num(mat, nan=np.nanmean(mat))
    H, W = X.shape

    # bands are fixed by design: even_k=2, odd_k=1
    even_k = 2
    odd_k = 1

    if square_size % 2 == 0:
        square_size += 1
    half = square_size // 2

    centers = np.arange(half, min(H, W) - half, 1, dtype=int)

    boxes = []
    for mid in centers:
        even_vals = band_vals_in_square(X, even_k, mid, half)
        odd_vals  = band_vals_in_square(X, odd_k,  mid, half)

        even_keep = trim_top(even_vals, trim_even_top)
        odd_keep  = trim_bottom(odd_vals, trim_odd_bottom)

        if even_keep.size > 0 and odd_keep.size > 0:
            max_even = np.nanmax(even_keep)
            min_odd  = np.nanmin(odd_keep)
            if np.isfinite(max_even) and np.isfinite(min_odd) and (max_even < min_odd):
                r1 = max(0, mid - half)
                r2 = min(H - 1, mid + half)
                c1 = r1
                c2 = r2
                boxes.append((int(r1), int(c1), int(r2), int(c2)))

    merged_boxes = merge_overlapping_or_touching(boxes)
    return merged_boxes


# === New: extract band values inside an arbitrary rectangle ===================

def band_vals_in_rect(mat, k, r1, c1, r2, c2):
    """
    Values along band offset +k limited to rectangle [r1..r2] x [c1..c2].
    Points are (r, r+k) that fall inside the rectangle.
    """
    H, W = mat.shape
    r1 = max(0, min(r1, H-1))
    r2 = max(0, min(r2, H-1))
    c1 = max(0, min(c1, W-1))
    c2 = max(0, min(c2, W-1))
    if r2 < r1 or c2 < c1:
        return np.array([], dtype=float)

    # rows such that r in [r1,r2] and r+k in [c1,c2]
    r_start = max(r1, c1 - k)
    r_end   = min(r2, c2 - k)
    if r_end < r_start:
        return np.array([], dtype=float)
    rows = np.arange(r_start, r_end + 1, dtype=int)
    cols = rows + k
    vals = mat[rows, cols]
    return vals.astype(float)


# ===============================
# Plotting
# ===============================

def plot_heatmap(dist_matrix, nonoverlapping_df, diagonals_df, output_png, invert_colors=False,
                 min_diag_len=4, diag_offset=0, diag_width=1.5, diag_color="black",
                 dist_min=-1, dist_max=-1, dpi=300, checker_zones=None):
    """Plot heatmap with diagonals and optional checkerboard zones."""
    cmap_name = "coolwarm_r" if not invert_colors else "coolwarm"
    cmap = plt.colormaps[cmap_name]

    fig, ax = plt.subplots(figsize=(10, 10))
    nrows, ncols = dist_matrix.shape
    fig_width, fig_height = fig.get_size_inches()
    min_dpi = max(ncols / fig_width, nrows / fig_height)
    dpi = int(np.ceil(max(dpi, min_dpi * 1.1)))

    if dist_min == -1:
        im = ax.imshow(dist_matrix, origin="upper", cmap=cmap, aspect="equal")
    else:
        im = ax.imshow(dist_matrix, origin="upper", cmap=cmap, aspect="equal",
                       vmin=dist_min, vmax=dist_max)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label("Distance", rotation=270, labelpad=15)
    cbar.ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    cbar.ax.tick_params(labelsize=8, length=3)

    # All diagonals (dashed) in upper triangle
    if diagonals_df is not None and not diagonals_df.empty:
        for _, row in diagonals_df.iterrows():
            x = np.arange(int(row["start_col"]), int(row["end_col"]) + 1) + diag_offset
            y = np.arange(int(row["start_row"]), int(row["end_row"]) + 1)
            mask = x >= y
            if np.any(mask):
                ax.plot(x[mask], y[mask], linestyle="--", color=diag_color,
                        lw=max(0.8, diag_width * 0.6), alpha=0.6)

    # Non-overlapping (solid) in upper triangle
    if nonoverlapping_df is not None and not nonoverlapping_df.empty:
        for _, row in nonoverlapping_df.iterrows():
            if row["length"] >= min_diag_len:
                x = np.arange(int(row["start_col"]), int(row["end_col"]) + 1) + diag_offset
                y = np.arange(int(row["start_row"]), int(row["end_row"]) + 1)
                mask = x >= y
                if np.any(mask):
                    ax.plot(x[mask], y[mask], color=diag_color, lw=diag_width)
                    mid_x = int(row["mid_col"]) + diag_offset
                    mid_y = int(row["mid_row"])
                    if mid_x >= mid_y:
                        ax.text(mid_x, mid_y, str(int(row["length"])),
                                color="black", fontsize=8, ha="center", va="center",
                                bbox=dict(facecolor="white", alpha=0.6, edgecolor="none", pad=1))

    # Checkerboard zones (rectangles)
    if checker_zones:
        for (r1, c1, r2, c2) in checker_zones:
            rect = plt.Rectangle((c1, r1), c2 - c1, r2 - r1,
                                 edgecolor="blue", facecolor="none", lw=2)
            ax.add_patch(rect)

    ax.set_xlabel("Columns")
    ax.set_ylabel("Rows")
    ax.set_title("Distance matrix with diagonals and checkerboard zones")

    plt.tight_layout()
    plt.savefig(output_png, dpi=dpi)
    plt.close()


# ===============================
# Worker (one file)
# ===============================

def _df_to_mask(df: pd.DataFrame, shape):
    """Build a boolean mask from a diagonals DataFrame."""
    M = np.zeros(shape, dtype=bool)
    for _, row in df.iterrows():
        r0, r1 = int(row["start_row"]), int(row["end_row"])
        c0, c1 = int(row["start_col"]), int(row["end_col"])
        rs = np.arange(r0, r1 + 1, dtype=int)
        cs = np.arange(c0, c1 + 1, dtype=int)
        M[rs, cs] = True
    return M


def _segments_to_df(segs, dist_matrix):
    """Turn list of Nx2 index arrays into a diagonals DataFrame."""
    recs = []
    for seg in segs:
        ii = seg[:, 0].astype(int)
        jj = seg[:, 1].astype(int)
        diag_id = int(jj[0] - ii[0])
        vals = dist_matrix[ii, jj]
        recs.append({
            "diag_id": diag_id,
            "start_row": int(ii.min()),
            "end_row": int(ii.max()),
            "start_col": int(jj.min()),
            "end_col": int(jj.max()),
            "length": int(len(ii)),
            "mid_row": int(np.round((ii.min() + ii.max()) / 2)),
            "mid_col": int(np.round((jj.min() + jj.max()) / 2)),
            "mean_distance": float(np.nanmean(vals))
        })
    return pd.DataFrame.from_records(recs) if recs else pd.DataFrame(columns=[
        "diag_id","start_row","end_row","start_col","end_col","length","mid_row","mid_col","mean_distance"
    ])


# === Helpers for combined & individual TSV outputs ============================

def _append_tsv(out_dir: str, filename: str, df: pd.DataFrame, block_id: str):
    """Append DataFrame to TSV in out_dir/filename, adding a leading 'block_id' column."""
    if df is None or (hasattr(df, "empty") and df.empty):
        return
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, filename)
    tmp = df.copy().reset_index(drop=True)
    # ensure block_id is first column
    tmp.insert(0, "block_id", block_id)
    header = not os.path.exists(out_path)
    tmp.to_csv(out_path, sep="\t", index=False, header=header, mode="a", na_rep="NA")


def _ensure_combined_header(out_dir: str, filename: str, columns: list):
    """Ensure the combined TSV exists with a header (no rows)."""
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, filename)
    if not os.path.exists(out_path):
        pd.DataFrame(columns=columns).to_csv(out_path, sep="\t", index=False, header=True, mode="w", na_rep="NA")


def _format_diags_for_combined(df: pd.DataFrame, nmonomers: int) -> pd.DataFrame:
    """
    Reformat diagonals for combined TSV:
      - drop: diagline_mean, region_median_full, region_mean_full, region_mean_excl,
              mean_distance, mid_row, mid_col
      - rename: diag_id -> order; region_median_excl -> median_dist_intra;
                diagline_median -> median_dist_inter
      - new diag_id string: start_row_end_row_start_col_end_col (keeps coords)
      - append monomers_count
    If empty: return a single NA row (best-effort) so each block appears at least once.
    """
    keep_cols = [
        "start_row", "end_row", "start_col", "end_col",
        "length", "region_median_excl", "diagline_median", "diag_id"
    ]

    if df is None or df.empty:
        out = pd.DataFrame([{
            "start_row": np.nan, "end_row": np.nan, "start_col": np.nan, "end_col": np.nan,
            "length": np.nan, "region_median_excl": np.nan, "diagline_median": np.nan,
            "diag_id": np.nan  # becomes "order"
        }])
    else:
        # ensure needed stats columns exist
        for c in ["region_median_excl", "diagline_median"]:
            if c not in df.columns:
                df[c] = np.nan
        out = df.copy()[[c for c in keep_cols if c in df.columns]]

    out = out.rename(columns={
        "region_median_excl": "median_dist_intra",
        "diagline_median": "median_dist_inter",
        "diag_id": "order",
    })

    def _mk_id(row):
        if pd.isna(row.get("start_row")) or pd.isna(row.get("end_row")) \
           or pd.isna(row.get("start_col")) or pd.isna(row.get("end_col")):
            return np.nan
        return f"{int(row['start_row'])}_{int(row['end_row'])}_{int(row['start_col'])}_{int(row['end_col'])}"
    out["diag_id"] = out.apply(_mk_id, axis=1)

    out["monomers_count"] = nmonomers

    desired = [
        "diag_id", "start_row", "end_row", "start_col", "end_col",
        "order", "length", "median_dist_intra", "median_dist_inter", "monomers_count"
    ]
    out = out[[c for c in desired if c in out.columns]]
    return out


def _format_checker_for_combined(df: pd.DataFrame, nmonomers: int) -> pd.DataFrame:
    """
    Format checkerboard zones for combined TSV with columns:
      block_id, checkerboard_id (r1_r2), start (r1), end (r2),
      median_dist_intra (band1_median), median_dist_inter (band2_median), monomers_count
    Drop c1, c2, band1_mean, band2_mean.
    If empty: return empty DataFrame (header ensured separately).
    """
    final_cols = ["checkerboard_id", "start", "end",
                  "median_dist_intra", "median_dist_inter", "monomers_count"]

    if df is None or df.empty:
        return pd.DataFrame(columns=final_cols)

    df2 = df.copy()
    rename_map = {"r1": "start", "r2": "end",
                  "band1_median": "median_dist_intra",
                  "band2_median": "median_dist_inter"}
    for old, new in rename_map.items():
        if old in df2.columns:
            df2 = df2.rename(columns={old: new})
        else:
            df2[new] = pd.Series(dtype=float)

    def _mk_id(row):
        if pd.isna(row.get("start")) or pd.isna(row.get("end")):
            return np.nan
        return f"{int(row['start'])}_{int(row['end'])}"
    df2["checkerboard_id"] = df2.apply(_mk_id, axis=1)

    for col in ["c1", "c2", "band1_mean", "band2_mean"]:
        if col in df2.columns:
            df2 = df2.drop(columns=[col])

    df2["monomers_count"] = nmonomers
    df2 = df2[final_cols]
    return df2


def _ensure_header_file(path: str, columns: list):
    """Create an empty TSV with header if it doesn't exist (or ensure header exists)."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    if not os.path.exists(path):
        pd.DataFrame(columns=columns).to_csv(path, sep="\t", index=False, header=True, mode="w", na_rep="NA")


def _format_diags_for_individual(df: pd.DataFrame, nmonomers: int) -> pd.DataFrame:
    """
    Same schema as combined diagonals but WITHOUT block_id.
    Columns (order):
      diag_id, start_row, end_row, start_col, end_col, order, length,
      median_dist_intra, median_dist_inter, monomers_count
    """
    cols = ["diag_id", "start_row", "end_row", "start_col", "end_col",
            "order", "length", "median_dist_intra", "median_dist_inter", "monomers_count"]
    if df is None or df.empty:
        return pd.DataFrame(columns=cols)

    keep_cols = [
        "start_row", "end_row", "start_col", "end_col",
        "length", "region_median_excl", "diagline_median", "diag_id"
    ]

    for c in ["region_median_excl", "diagline_median"]:
        if c not in df.columns:
            df[c] = np.nan

    out = df.copy()[[c for c in keep_cols if c in df.columns]].rename(columns={
        "region_median_excl": "median_dist_intra",
        "diagline_median": "median_dist_inter",
        "diag_id": "order",
    })

    def _mk_id(row):
        if pd.isna(row.get("start_row")) or pd.isna(row.get("end_row")) \
           or pd.isna(row.get("start_col")) or pd.isna(row.get("end_col")):
            return np.nan
        return f"{int(row['start_row'])}_{int(row['end_row'])}_{int(row['start_col'])}_{int(row['end_col'])}"
    out["diag_id"] = out.apply(_mk_id, axis=1)

    out["monomers_count"] = nmonomers
    out = out[cols]
    return out


def _format_checker_for_individual(df: pd.DataFrame, nmonomers: int) -> pd.DataFrame:
    """
    Same schema as combined checkerboard but WITHOUT block_id.
    Columns (order):
      checkerboard_id, start, end, median_dist_intra, median_dist_inter, monomers_count
    """
    cols = ["checkerboard_id", "start", "end", "median_dist_intra", "median_dist_inter", "monomers_count"]
    if df is None or df.empty:
        return pd.DataFrame(columns=cols)

    df2 = df.copy()
    rename_map = {"r1": "start", "r2": "end",
                  "band1_median": "median_dist_intra",
                  "band2_median": "median_dist_inter"}
    for old, new in rename_map.items():
        if old in df2.columns:
            df2 = df2.rename(columns={old: new})
        else:
            df2[new] = pd.Series(dtype=float)

    def _mk_id(row):
        if pd.isna(row.get("start")) or pd.isna(row.get("end")):
            return np.nan
        return f"{int(row['start'])}_{int(row['end'])}"
    df2["checkerboard_id"] = df2.apply(_mk_id, axis=1)

    for col in ["c1", "c2", "band1_mean", "band2_mean"]:
        if col in df2.columns:
            df2 = df2.drop(columns=[col])

    df2["monomers_count"] = nmonomers
    return df2[cols]


def _process_one_file(path, cfg):
    """Process a single matrix file; returns (path, ok, error_message_or_None)."""
    try:
        dist_matrix = load_matrix(path, cfg["input_type"])
        stem = os.path.splitext(os.path.basename(path))[0]
        nmonomers = int(dist_matrix.shape[0])

        # Determine prefix (batch-aware)
        if cfg["batch_mode"]:
            if cfg["output_prefix"]:
                prefix = f"{cfg['output_prefix']}_{stem}"
            else:
                prefix = stem
        else:
            prefix = cfg["output_prefix"] or stem

        # Output paths
        out_dir = cfg["out_dir"]
        diag_csv = os.path.join(out_dir, f"{prefix}_diagonals.tsv")
        nonov_csv = os.path.join(out_dir, f"{prefix}_nonoverlapping_diagonals.tsv")
        heat_png = os.path.join(out_dir, f"{prefix}_heatmap.png")
        checker_csv = os.path.join(out_dir, f"{prefix}_checkerboard_zones.tsv")

        # Ensure per-file TSVs exist with header even if no detections
        _ensure_header_file(diag_csv, ["diag_id","start_row","end_row","start_col","end_col","order","length","median_dist_intra","median_dist_inter","monomers_count"])
        _ensure_header_file(nonov_csv, ["diag_id","start_row","end_row","start_col","end_col","order","length","median_dist_intra","median_dist_inter","monomers_count"])
        _ensure_header_file(checker_csv, ["checkerboard_id","start","end","median_dist_intra","median_dist_inter","monomers_count"])

        diagonals_df, nonoverlapping_df = pd.DataFrame(), pd.DataFrame()
        checker_zones = []
        df_boxes = None

        if cfg["detect_diagonals"]:
            # 1) Detection (+ merge if requested)
            diagonals_df = extract_diagonals(
                dist_matrix,
                conv_size=cfg["conv_size"],
                local_size=cfg["local_size"],
                sd_factor=cfg["sd_factor"],
                min_length=cfg["min_diag_len"],
                merge=cfg["merge_diagonals"],
                merge_gap=cfg["merge_gap"],
                merge_mode=cfg["merge_mode"]
            )

            # 2) Optional peeling (upper triangle)
            if cfg["peel"] and not diagonals_df.empty:
                mask0 = _df_to_mask(diagonals_df, dist_matrix.shape)
                side = "upper"
                _, segs_after = peel_mask_by_contrast(
                    dist_matrix, mask0,
                    side=side,
                    k_offsets=cfg["peel_k_offsets"],
                    min_len=cfg["min_diag_len"],
                    max_peel_per_side=cfg["peel_max_per_side"],
                    combine=cfg["peel_combine"],
                )
                diagonals_df = _segments_to_df(segs_after, dist_matrix)

            # 3) Stats per-file (we only need medians for final TSVs)
            if not diagonals_df.empty:
                reg_cols = [
                    "region_mean_excl", "region_median_excl",
                    "region_mean_full", "region_median_full",
                    "diagline_mean", "diagline_median"
                ]
                region_summaries = diagonals_df.apply(
                    lambda row: pd.Series(compute_region_summaries(dist_matrix, row), index=reg_cols),
                    axis=1
                )
                diagonals_df = pd.concat([diagonals_df, region_summaries], axis=1)

            # 4) Non-overlapping + per-file TSVs
            nonoverlapping_df = select_non_overlapping(diagonals_df, cfg["diag_select"]).round(2)
            _format_diags_for_individual(diagonals_df.round(2), nmonomers).to_csv(diag_csv, index=False, sep="\t", na_rep="NA")
            _format_diags_for_individual(nonoverlapping_df, nmonomers).to_csv(nonov_csv, index=False, sep="\t", na_rep="NA")

            # 5) Combined TSVs for diagonals
            if cfg.get("combined_outputs"):
                diag_for_tsv  = _format_diags_for_combined(diagonals_df.round(2), nmonomers)
                nonov_for_tsv = _format_diags_for_combined(nonoverlapping_df, nmonomers)
                _append_tsv(cfg["combined_outputs"], "combined_diagonals.tsv", diag_for_tsv, stem)
                _append_tsv(cfg["combined_outputs"], "combined_nonoverlapping_diagonals.tsv", nonov_for_tsv, stem)

        # 6) Checkerboard detection (parity-trim)
        if cfg["detect_checkerboard"]:
            checker_zones = detect_checkerboard_parity_trimmed_sliding(
                dist_matrix,
                square_size=cfg["checker_square"],
                trim_even_top=cfg["checker_trim_even_top"],
                trim_odd_bottom=cfg["checker_trim_odd_bottom"],
            )

            # Build dataframe with stats per merged square
            cols = ["r1", "c1", "r2", "c2",
                    "band1_mean", "band1_median",
                    "band2_mean", "band2_median"]
            rows = []
            if checker_zones:
                for (r1, c1, r2, c2) in checker_zones:
                    vals_b1 = band_vals_in_rect(dist_matrix, 1, r1, c1, r2, c2)
                    vals_b2 = band_vals_in_rect(dist_matrix, 2, r1, c1, r2, c2)

                    b1_mean = float(np.nanmean(vals_b1)) if vals_b1.size else np.nan
                    b1_med  = float(np.nanmedian(vals_b1)) if vals_b1.size else np.nan
                    b2_mean = float(np.nanmean(vals_b2)) if vals_b2.size else np.nan
                    b2_med  = float(np.nanmedian(vals_b2)) if vals_b2.size else np.nan

                    rows.append([int(r1), int(c1), int(r2), int(c2),
                                 b1_mean, b1_med, b2_mean, b2_med])

                df_boxes = pd.DataFrame(rows, columns=cols).round(2)
            else:
                df_boxes = pd.DataFrame(columns=cols)

            # Per-file TSV
            _format_checker_for_individual(df_boxes, nmonomers).to_csv(checker_csv, index=False, sep="\t", na_rep="NA")

            # Combined TSV: ensure header always exists, append rows only if any
            if cfg.get("combined_outputs"):
                _ensure_combined_header(cfg["combined_outputs"], "combined_checkerboard_zones.tsv",
                                        ["block_id", "checkerboard_id", "start", "end",
                                         "median_dist_intra", "median_dist_inter", "monomers_count"])
                cb_for_tsv = _format_checker_for_combined(df_boxes, nmonomers)
                if cb_for_tsv is not None and not cb_for_tsv.empty:
                    _append_tsv(cfg["combined_outputs"], "combined_checkerboard_zones.tsv", cb_for_tsv, stem)

        # 7) Heatmap with overlays
        plot_heatmap(
            dist_matrix, nonoverlapping_df, diagonals_df, heat_png,
            invert_colors=cfg["invert_colors"],
            min_diag_len=cfg["min_diag_len"],
            diag_offset=cfg["diag_offset"],
            diag_color=cfg["diag_color"],
            diag_width=cfg["diag_width"],
            dist_min=cfg["vmin"],
            dist_max=cfg["vmax"],
            dpi=cfg["dpi"],
            checker_zones=checker_zones
        )
        return (path, True, None)
    except Exception as e:
        return (path, False, str(e))


# ===============================
# Main
# ===============================

def _parse_k_offsets(s: str):
    """Parse '1,2' -> [1,2] ; ignore spaces."""
    parts = [p.strip() for p in s.split(",") if p.strip() != ""]
    out = []
    for p in parts:
        out.append(int(p))
    return out if out else [1, 2]


def main():
    parser = argparse.ArgumentParser(
        description="Detect diagonals and checkerboard in distance matrices",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Input: single file OR directory
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--input", help="Single input distance matrix file (.txt or .h5)")
    group.add_argument("--input_dir", help="Directory containing multiple matrices")

    parser.add_argument("--input_type", choices=["text", "h5"], default="h5", help="Input file type")
    parser.add_argument("--output_prefix", type=str, default=None,
                        help="Optional prefix for output files. In batch mode: <prefix>_<file_stem>_*; otherwise: <prefix>_*.")
    parser.add_argument("--output_dir", type=str, default=None,
                        help="Optional output directory. Defaults to input file's directory or to --input_dir in batch mode.")
    parser.add_argument("--workers", type=int, default=os.cpu_count(),
                        help="Number of parallel workers (processes) to use in batch mode")
    parser.add_argument("--invert_colors", action="store_true", help="Invert heatmap colors")
    parser.add_argument("--vmin", type=float, default=-1, help="Minimum value for heatmap color scale")
    parser.add_argument("--vmax", type=float, default=-1, help="Maximum value for heatmap color scale")
    parser.add_argument("--dpi", type=int, default=300, help="DPI for heatmap (auto-adjusts if too low)")

    # Diagonals
    parser.add_argument("--detect_diagonals", action="store_true", help="Enable diagonal detection")
    parser.add_argument("--conv_size", type=int, default=5, help="Diagonal convolution kernel size (odd)")
    parser.add_argument("--local_size", type=int, default=7, help="Local window size for robust stats")
    parser.add_argument("--sd_factor", type=float, default=2.0, help="Robust SD factor for detection")
    parser.add_argument("--min_diag_len", type=int, default=5, help="Minimum diagonal length to annotate")
    parser.add_argument("--merge_diagonals", action="store_true", help="Enable merging of diagonal fragments")
    parser.add_argument("--merge_gap", type=int, default=5, help="Gap for merging diagonals")
    parser.add_argument("--merge_mode", choices=["fixed", "relative"], default="fixed", help="Merge mode")
    parser.add_argument("--diag_select", choices=["repet", "length", "distance"], default="length",
                        help="Criterion for selecting non-overlapping diagonals: 'repet'=length/diag_id, 'length'=longest, 'distance'=smallest mean distance")
    parser.add_argument("--diag_offset", type=int, default=0, help="Horizontal offset when plotting diagonals")
    parser.add_argument("--diag_color", type=str, default="black", help="Color of diagonal lines")
    parser.add_argument("--diag_width", type=float, default=1.5, help="Line width of diagonal lines")

    # Peeling
    parser.add_argument("--peel", action="store_true", help="Enable end-peeling of diagonal segments by parallel-background contrast")
    parser.add_argument("--peel_k_offsets", type=str, default="1,2", help="Parallel diagonal offsets to build background (e.g. '1,2')")
    parser.add_argument("--peel_combine", type=str, default="p10", help="Combine rule for background: median|mean|pXX|qXX (default: p10)")
    parser.add_argument("--peel_max_per_side", type=int, default=None, help="Max points to peel at head/tail (None = unlimited)")

    # Checkerboard (parity + trim)
    parser.add_argument("--detect_checkerboard", action="store_true", help="Enable checkerboard detection by parity-trim in sliding squares")
    parser.add_argument("--checker_square", type=int, default=30, help="Sliding square side (objects). If even, it is forced to odd by +1.")
    parser.add_argument("--checker_trim_even_top", type=float, default=0.20, help="Trim top fraction on even band inside square (e.g., 0.20)")
    parser.add_argument("--checker_trim_odd_bottom", type=float, default=0.20, help="Trim bottom fraction on odd band inside square (e.g., 0.20)")

    # Combined aggregated outputs (TSV append)
    parser.add_argument("--combined_outputs", type=str, default=None,
                        help="Directory where to append aggregated TSVs: combined_diagonals.tsv, combined_nonoverlapping_diagonals.tsv, combined_checkerboard_zones.tsv; each row prefixed with block_id=<file stem>")

    args = parser.parse_args()

    # Prepare list of files to process
    if args.input:
        files = [args.input]
        default_out_dir = os.path.dirname(os.path.abspath(args.input)) or "."
        batch_mode = False
    else:
        if not os.path.isdir(args.input_dir):
            raise FileNotFoundError(f"Input directory not found: {args.input_dir}")
        ext = ".h5" if args.input_type == "h5" else ".txt"
        files = sorted(
            os.path.join(args.input_dir, f) for f in os.listdir(args.input_dir)
            if f.lower().endswith(ext)
        )
        if not files:
            raise FileNotFoundError(f"No '{ext}' files found in: {args.input_dir}")
        default_out_dir = os.path.abspath(args.input_dir)
        batch_mode = True

    out_dir = os.path.abspath(args.output_dir) if args.output_dir else default_out_dir
    os.makedirs(out_dir, exist_ok=True)

    cfg = {
        "input_type": args.input_type,
        "output_prefix": args.output_prefix,
        "out_dir": out_dir,
        "batch_mode": batch_mode,
        "invert_colors": args.invert_colors,
        "vmin": args.vmin,
        "vmax": args.vmax,
        "dpi": args.dpi,

        "detect_diagonals": args.detect_diagonals,
        "conv_size": args.conv_size,
        "local_size": args.local_size,
        "sd_factor": args.sd_factor,
        "min_diag_len": args.min_diag_len,
        "merge_diagonals": args.merge_diagonals,
        "merge_gap": args.merge_gap,
        "merge_mode": args.merge_mode,
        "diag_select": args.diag_select,
        "diag_offset": args.diag_offset,
        "diag_color": args.diag_color,
        "diag_width": args.diag_width,

        "peel": args.peel,
        "peel_k_offsets": _parse_k_offsets(args.peel_k_offsets),
        "peel_combine": args.peel_combine,
        "peel_max_per_side": args.peel_max_per_side,

        "detect_checkerboard": args.detect_checkerboard,
        "checker_square": args.checker_square,
        "checker_trim_even_top": args.checker_trim_even_top,
        "checker_trim_odd_bottom": args.checker_trim_odd_bottom,

        # combined outputs dir (or None)
        "combined_outputs": (os.path.abspath(args.combined_outputs)
                             if args.combined_outputs else None),
    }

    # Process (parallel in batch mode if workers > 1)
    if batch_mode and args.workers and args.workers > 1 and len(files) > 1:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futures = {ex.submit(_process_one_file, path, cfg): path for path in files}
            for fut in as_completed(futures):
                path = futures[fut]
                ok = False
                err = None
                try:
                    _, ok, err = fut.result()
                except Exception as e:
                    ok, err = False, str(e)
                if not ok:
                    print(f"[ERROR] {os.path.basename(path)}: {err}")
                else:
                    print(f"[OK]    {os.path.basename(path)}")
    else:
        for path in files:
            _, ok, err = _process_one_file(path, cfg)
            if not ok:
                print(f"[ERROR] {os.path.basename(path)}: {err}")
            else:
                print(f"[OK]    {os.path.basename(path)}")


if __name__ == "__main__":
    main()
