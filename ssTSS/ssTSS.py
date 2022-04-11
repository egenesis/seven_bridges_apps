import HTSeq
import itertools
import pandas as pd
import argparse
import os

## Parse Arguments

parser = argparse.ArgumentParser(description='Extract TSSs')
parser.add_argument('-b', '--bam', type=str, required=True, help='Input BAM')
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

# Remove duplicates.
TSSs.sort()
TSSs = list(x for x, _ in itertools.groupby(TSSs))

# Turn into Pandas DataFrame and save as csv.
df = pd.DataFrame(TSSs, columns=['chrm', 'strand', 'tss', 'cell_barcode', 'umi'])
df.to_csv(os.path.splitext(args.bam)[0] + '.csv.gz', index=False, compression='gzip')