#!/usr/bin/env python3

import HTSeq
import itertools
import pandas as pd
import argparse
import os
import numpy as np
import scipy.sparse as sparse
import scipy.io as sio
from datetime import datetime
import gc

#####################
## Parse Arguments ##
#####################

parser = argparse.ArgumentParser(description='Extract TSSs')
parser.add_argument('-b', '--bam', type=str, required=True, help='Input BAM')
parser.add_argument('-n', '--ncells', type=int, required=True, help='Number of cells required to keep a TSS')
args = parser.parse_args()


## Test values.

#class MakeArgs(object):
#    def __init__(self, bam, ncells):
#        self.bam = bam
#        self.ncells = ncells

#args = MakeArgs('', 10)

###################
## Retrieve TSSs ##
###################

TSSs = []

print(f"{datetime.now()} - Reading {args.bam}")
with HTSeq.BAM_Reader(args.bam) as f:
    for (read1, read2) in HTSeq.pair_SAM_alignments(f):
        # Store some read and alignment information.
        try:
            cell_barcode = read1.optional_field('CB')
            umi = read1.optional_field('UB')
        except AttributeError:
            continue
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
        # Add TSS to list.
        TSS = (f"{chrm}:{tss}:{strand}", cell_barcode, umi)
        TSSs.append(TSS)

gc.collect()

##################
## Process TSSs ##
##################

# Remove duplicates within cells.
print(f"{datetime.now()} - Removing duplicates")
TSSs = list(x[0:2] for x, _ in itertools.groupby(TSSs.sort()))

# Aggregate overlapping counts.
print(f"{datetime.now()} - Aggregating overlapping counts")
TSSs = [(x[0], x[1], y) for x, y in Counter(TSSs).items()]

# Convert to DataFrame.
print(f"{datetime.now()} - Converting to a pandas DataFrame")
TSSs = pd.DataFrame(TSSs, columns=['range', 'cell_barcode', 'score'])

# Convert the score column to integers.
print(f"{datetime.now()} - Converting score to integers")
TSSs['score'] = TSSs['score'].astype(int)

# Find number of cells that a TSS is present in.
print(f"{datetime.now()} - Finding number of cells that a TSS is present in")
TSSs['n_cells'] = TSSs.groupby('range')['range'].transform('size')

# Filter out TSSs that are not present in enough cells.
print(f"{datetime.now()} - Filtering out TSSs that are not present in at least {args.ncells} cells")
TSSs = TSSs[TSSs['n_cells'] >= args.ncells]

############################
## Convert to Wide Format ##
############################

# Get unique barcodes and ranges.
barcodes = list(set(TSSs['cell_barcode'])).sort()
ranges = list(set(TSSs['range']))

# Create the counts matrix.
print(f"{datetime.now()} - Creating the counts matrix")
def fill_mat(row):
    tss_mat[ranges.index(row['range']), barcodes.index(row['cell_barcode'])] = row['score']

tss_mat = np.zeros((len(ranges), len(barcodes)))

print(f"{datetime.now()} - Filling the counts matrix")
TSSs.apply(lambda row: fill_mat(row), axis=1)

# Convert the matrix to a sparse matrix.
tss_mat = sparse.csr_matrix(tss_mat)

#####################
## Export the TSSs ##
#####################

# Create output directory.
output_dir = os.path.splitext(os.path.basename(args.bam))[0]
os.mkdir(output_dir)

# Save the barcodes.
print(f"{datetime.now()} - Saving the barcodes")
barcodes = pd.DataFrame(TSSs.columns)
barcodes.to_csv(f"{output_dir}/barcodes.tsv", index=False, header=False, sep="\t")

# Save the features.
print(f"{datetime.now()} - Saving the features")
features = pd.DataFrame(TSSs.index)
features.to_csv(f"{output_dir}/features.tsv", index=False, header=False, sep="\t")

# Save the matrix.
print(f"{datetime.now()} - Saving the matrix")
sio.mmwrite(f"{output_dir}/counts.mtx", tss_mat)
