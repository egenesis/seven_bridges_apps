#!/usr/bin/env python

import HTSeq
import collections
import argparse
import os
from datetime import datetime
import gc
import dask.dataframe as dd

#####################
## Parse Arguments ##
#####################

parser = argparse.ArgumentParser(description='Extract TSSs')
parser.add_argument('-b', '--bam', type=str, required=True, help='Input BAM')
parser.add_argument('-m', '--ncells', type=int, required=True, help='Minimum number of cells a TSS must be observed in')
parser.add_argument('-t', '--tcount', type=int, required=True, help='at least tcells cells must have a score of tcount')
parser.add_argument('-c', '--tcells', type=int, required=True, help='at least tcells cells must have a score of tcount')
parser.add_argument('-f', '--nfeatures', type=int, required=True, help='Number of features required to keep a cell')
args = parser.parse_args()


## Test values.

#class MakeArgs(object):
#    def __init__(self, bam, ncells, nfeatures, tcount, tcells):
#        self.bam = bam
#        self.ncells = ncells
#        self.nfeatures = nfeatures
#        self.tcount = tcount
#        self.tcells = tcells

#args = MakeArgs('', 10, 100, 3, 3)

###################
## Retrieve TSSs ##
###################

print(f"{datetime.now()} - Reading {args.bam}")

TSSs = collections.deque()
i = 0
with open("tempTSSs.csv", "w") as ftss:
    ftss.write("range,barcode,umi\n")
    with HTSeq.BAM_Reader(args.bam) as fbam:
        for (read1, read2) in HTSeq.pair_SAM_alignments(fbam):
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
            TSSs.append((f"{chrm}:{tss}:{strand},{cell_barcode},{umi}\n"))
            # Print number of reads processed for every 1E6 reads.
            i += 1
            if i % 1E6 == 0:
                print(f"{datetime.now()} - processed {i} valid read pairs")
            # For every 1E7 reads, save to file to keep memory usage low.
            if i % 1E7 == 0:
                print(f"{datetime.now()} - writing {len(TSSs)} TSSs to tempTSSs.csv")
                ftss.write("".join(TSSs))
                print(f"{datetime.now()} - finished writing {len(TSSs)} TSSs to tempTSSs.csv")
                TSSs = collections.deque()
    # Write any last TSSs to file.
    if i % 1E7 != 0:
        print(f"{datetime.now()} - writing {len(TSSs)} TSSs to tempTSSs.csv")
        ftss.write("".join(TSSs))
        print(f"{datetime.now()} - finished writing {len(TSSs)} TSSs to tempTSSs.csv")
    print(f"{datetime.now()} - wrote {i} total TSSs to tempTSSs.csv")

del TSSs
gc.collect()

##################
## Process TSSs ##
##################

# Create dask DataFrame.
print(f"{datetime.now()} - Creating the dask DataFrame")

df = dd.read_csv("tempTSSs.csv")

# Remove duplicates.
print(f"{datetime.now()} - Removing duplicates")

df = df.drop_duplicates().drop('umi', axis=1)

# Aggregate overlapping counts.
print(f"{datetime.now()} - Aggregating overlapping counts")

df = df.groupby(['range', 'barcode']).size().reset_index().rename(columns={0: 'score'})

# Remove ranges that are not present in at least ncells.
print(f"{datetime.now()} - Finding number of cells that a TSS is present in")

df['ncells'] = df.groupby(['range'])['barcode'].transform('nunique', meta=('barcode', 'int'))
df = df[df['ncells'] >= args.ncells].drop('ncells', axis=1)

# Remove remaining ranges that don't have at least tcells cells with a count of tcount.
print(f"{datetime.now()} - Finding number of cells that a TSS has a count of {args.tcount}")

filtered = df[df['score'] >= args.tcount].groupby('range')['barcode'].count().reset_index().rename(columns={'barcode': 'tcount'})
filtered = filtered[filtered['tcount'] >= args.tcells].drop('tcount', axis=1)
df = df.merge(filtered, how='inner', on='range')

del filtered

# Remove cells with less than nfeatures features.
print(f"{datetime.now()} - Finding number of features per cell")

df['nfeatures'] = df.groupby('barcode')['range'].transform('nunique', meta=('range', 'int'))
df = df[df['nfeatures'] >= args.nfeatures].drop('nfeatures', axis=1)

# Compute the results of the above filtering.
print(f"{datetime.now()} - Computing the filtering results")

df = df.compute()

# Add the row and column indices.
print(f"{datetime.now()} - Adding the row and column indices")

rowids = df['range'].drop_duplicates().to_frame()
rowids['rowid'] = list(range(1, len(rowids) + 1))

colids = df['barcode'].drop_duplicates().to_frame()
colids['colid'] = list(range(1, len(colids) + 1))

df = df.merge(rowids, 'left', 'range')
df = df.merge(colids, 'left', 'barcode')
df = df.sort_values(['rowid', 'colid'])

#####################
## Export the TSSs ##
#####################

# Create output directory.
print(f"{datetime.now()} - Create the output directory")

output_dir = os.path.splitext(os.path.basename(args.bam))[0]
os.mkdir(output_dir)

# Save the features.
print(f"{datetime.now()} - Saving the features")

rowids['range'].to_csv(f"{output_dir}/features.tsv", index=False, header=False)

# Save the barcodes.
print(f"{datetime.now()} - Saving the barcodes")

colids['barcode'].to_csv(f"{output_dir}/barcodes.tsv", index=False, header=False)

# Save the matrix.
print(f"{datetime.now()} - Saving the matrix")

with open(f"{output_dir}/matrix.mtx", "w") as mf:
    mf.write("%%MatrixMarket matrix coordinate integer general\n%\n")
    mf.write(f"{len(rowids)} {len(colids)} {len(df)}\n")

df[['rowid', 'colid', 'score']].to_csv(f"{output_dir}/matrix.mtx", index=False, header=False, sep=" ", mode="a")

print(f"{datetime.now()} - Finished!")