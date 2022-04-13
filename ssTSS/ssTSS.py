#!/usr/bin/env python3

import HTSeq
import itertools
import pandas as pd
import argparse
import os
import numpy as np
import scipy.sparse as sparse
import scipy.io as sio

## Parse Arguments

parser = argparse.ArgumentParser(description='Extract TSSs')
parser.add_argument('-b', '--bam', type=str, required=True, help='Input BAM')
parser.add_argument('-n', '--ncells', type=int, required=True, help='Number of cells required to keep a TSS')
args = parser.parse_args()

## Retrieve TSSs

TSSs = []

with HTSeq.BAM_Reader(args.bam) as f:
    for (read1, read2) in HTSeq.pair_SAM_alignments_with_buffer(f):
        # Check for some read pair issues.
        if not read1.aligned or not read2.aligned:
            continue
        if read1.not_primary_alignment or read2.not_primary_alignment:
            continue
        if read1.pcr_or_optical_duplicate or read2.pcr_or_optical_duplicate:
            continue
        if read1.supplementary or read2.supplementary:
            continue
        if not read1.proper_pair:
            continue
        # Store some read and alignment information.
        cell_barcode = read1.optional_field('CB')
        umi = read1.optional_field('UB')
        # Skip if no cell barcode or UMI present.
        if cell_barcode == "-" or umi == "-":
            continue
        # Get some genomic interval locations.
        chrm = read1.iv.chrom
        strand = read1.iv.strand
        # Get the TSS.
        if read1.pe_which == 'first':
            tss = read1.iv.start_d
        else:
            tss = read2.iv.start_d
        # Add TSS to list if not already present.
        TSS = (chrm, strand, tss, cell_barcode, umi)
        TSSs.append(TSS)

## Process TSSs

# Remove duplicates within cells.
TSSs.sort()
TSSs = list(x for x, _ in itertools.groupby(TSSs))

# Convert to DataFrame.
TSSs = [(f"{str(x[0])}:{str(x[2])}:{str(x[1])}", x[3]) for x in TSSs]
TSSs = pd.DataFrame(TSSs, columns=['range', 'cell_barcode'])

# Aggregate overlapping TSSs.
TSSs = TSSs.groupby(['range', 'cell_barcode']).size().reset_index(name="score")

# Convert the score column to integers.
TSSs['score'] = TSSs['score'].astype(int)

# Find number of cells that a TSS is present in.
TSSs['n_cells'] = TSSs.groupby('range')['range'].transform('size')

# Filter out TSSs that are not present in enough cells.
TSSs = TSSs[TSSs['n_cells'] >= args.ncells]

# Convert to wide format.
TSSs = TSSs.groupby(['range', 'cell_barcode'])['score'].max().unstack().replace(np.nan, 0)

## Export the TSSs

# Create output directory.
output_dir = os.path.splitext(os.path.basename(args.bam))[0]
os.mkdir(output_dir)

# Save the barcodes.
barcodes = pd.DataFrame(TSSs.columns)
barcodes.to_csv(f"{output_dir}/barcodes.tsv", index=False, header=False, sep="\t")

# Save the features.
features = pd.DataFrame(TSSs.index)
features.to_csv(f"{output_dir}/features.tsv", index=False, header=False, sep="\t")

# Save the matrix.
counts = sparse.csr_matrix(TSSs.to_numpy())
sio.mmwrite(f"{output_dir}/counts.mtx", counts)
