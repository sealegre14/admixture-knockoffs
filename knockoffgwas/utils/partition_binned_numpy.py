#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd

# -------------------------------
# Check arguments
# -------------------------------
if len(sys.argv) != 5:
    sys.exit(
        "Usage: python partition_binned_5cM_closest.py "
        "<map.file> <bim.file> <qc.file> <out.file>"
    )

map_file = sys.argv[1]
bim_file = sys.argv[2]
qc_file  = sys.argv[3]
out_file = sys.argv[4]

print(f"Map file: {map_file}")
print(f"BIM file: {bim_file}")
print(f"QC file: {qc_file}")
print(f"Output file: {out_file}\n")

# -------------------------------
# Load BIM (SNP positions)
# -------------------------------
bim = pd.read_csv(
    bim_file,
    sep="\t",
    header=None,
    names=["CHR", "SNP", "CM", "BP", "A1", "A2"]
)[["CHR", "SNP", "BP"]]

# -------------------------------
# Load QC SNP list
# -------------------------------
qc = pd.read_csv(
    qc_file,
    sep=r"\s+",
    header=None,
    names=["SNP"]
)

# Keep only QC-passed SNPs
# SA wrapping statment in if 
if qc.shape[0] > 0:
    bim = bim.merge(qc, on="SNP", how="inner")

# -------------------------------
# Load genetic map
# -------------------------------
gmap = pd.read_csv(
    map_file,
    sep="\t"
)[["pos", "chr", "cM"]]

# Ensure numeric + sorted
gmap["pos"] = gmap["pos"].astype(int)
gmap["cM"]  = gmap["cM"].astype(float)
gmap = gmap.sort_values("pos").reset_index(drop=True)

# -------------------------------
# Sort BIM by position
# -------------------------------
bim = bim.sort_values("BP").reset_index(drop=True)

# -------------------------------
# Match each SNP to closest map position
# -------------------------------
map_pos = gmap["pos"].values
map_cM  = gmap["cM"].values

snp_bp = bim["BP"].values

# Find insertion points
idx = np.searchsorted(map_pos, snp_bp)

# Clamp indices
idx_left  = np.clip(idx - 1, 0, len(map_pos) - 1)
idx_right = np.clip(idx,     0, len(map_pos) - 1)

# Compare distances
dist_left  = np.abs(snp_bp - map_pos[idx_left])
dist_right = np.abs(snp_bp - map_pos[idx_right])

use_right = dist_right < dist_left
closest_idx = np.where(use_right, idx_right, idx_left)

# Assign cM
bim["cM"] = map_cM[closest_idx]

# -------------------------------
# Bin into fixed 5 cM bins
# -------------------------------

# BIN_SIZE = 5.0
#bim["res_1"] = (bim["cM"] // BIN_SIZE).astype(int) + 1

# SA im gonna try out multiple binning sizes
BIN_SIZES = [1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 12.0]

# Create bins for each resolution
res_cols = {}
for i, bin_size in enumerate(BIN_SIZES):
    # res_1 = finest, res_k = coarsest (match partition.r convention)
    col_name = f"res_{len(BIN_SIZES) - i}"
    res_cols[col_name] = (bim["cM"] // bin_size).astype(int) + 1

# Add to dataframe
for col, values in res_cols.items():
    bim[col] = values



# -------------------------------
# Write output (knockoffgwas-compatible)
# -------------------------------

# SA I AM CHANGING THIS ENTIRE PART FOR NOW TO TRY WITH MULTIPLE RESOLUTIONS (cM SIZES)

# out = bim[["SNP", "res_1"]]

# out.to_csv(
#     out_file,
#     sep=" ",
#     index=False
# )

# print(f"Partitioning complete.")
# print(f"Number of SNPs: {len(out)}")
# print(f"Number of bins: {out['res_1'].nunique()}")
# print(f"Output written to: {out_file}")


'''new write output here'''

# Collect resolution columns in numeric order
res_cols = sorted(
    [c for c in bim.columns if c.startswith("res_")],
    key=lambda x: int(x.split("_")[1]),
    reverse=True  # <- coarsest first
)

out = bim[["SNP"] + res_cols]

# out.to_csv(out_file, sep=" ", index=False)
out.to_csv(out_file, sep=" ", index=False, header=True, float_format="%.0f")


print(f"Partitioning complete.")
print(f"Number of SNPs: {len(out)}")
for c in res_cols:
    print(f"{c}: {out[c].nunique()} bins")
print(f"Output written to: {out_file}")
