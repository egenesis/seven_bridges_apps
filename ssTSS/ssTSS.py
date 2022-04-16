#!/usr/bin/env spark-submit

import HTSeq
import itertools
import argparse
import os
from datetime import datetime
import gc
from pyspark.sql import SparkSession
from pyspark.sql import Window
import pyspark.sql.functions as funcs

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

#args = MakeArgs('hPBMC-E014_S1_L001_R1_001.fastq_Aligned.sortedByCoord.out.subsample.cleaned.bam', 20, 5, 3, 5)

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
TSSs.sort()
TSSs = list(x[0:2] for x, _ in itertools.groupby(TSSs))

# Create pyspark DataFrame.
spark = SparkSession.builder.appName('ssTSSs').getOrCreate()
df = spark.sparkContext.parallelize(TSSs)
df = df.toDF(['range', 'barcode'])

# Aggregate overlapping counts.
print(f"{datetime.now()} - Aggregating overlapping counts")
df = df.groupBy(['range', 'barcode']).count()

# Remove ranges that are not present in at least ncells.
print(f"{datetime.now()} - Finding number of cells that a TSS is present in")
w = Window.partitionBy('range')
df = df.select('range', 'barcode', 'count', funcs.count('range').over(w).alias('ncells'))
df = df[df.ncells >= args.ncells]
df = df.drop('ncells')

# Remove remaining ranges that don't have at least tcells cells with a count of tcount.
print(f"{datetime.now()} - Finding number of cells that a TSS has a count of {args.tcount}")
range_counts = df[df['count'] >= args.tcount].groupBy('range').count()
range_counts = range_counts[range_counts['count'] >= args.tcells].drop('count')
df = df.join(range_counts, 'range', how='inner')

# Remove cells with less than nfeatures features.
print(f"{datetime.now()} - Finding number of features per cell")
w = Window.partitionBy('barcode')
df = df.select('range', 'barcode', 'count', funcs.count('barcode').over(w).alias('nfeatures'))
df = df[df.nfeatures >= args.nfeatures].drop('nfeatures')

# Add the row and column indices.
rowids = df.select('range').distinct().sort('range').withColumn('rowid', funcs.monotonically_increasing_id())
rowids = rowids.withColumn('rowid', rowids.rowid + 1)

colids = df.select('barcode').distinct().sort('barcode').withColumn('colid', funcs.monotonically_increasing_id())
colids = colids.withColumn('colid', colids.colid + 1)

df = df.join(rowids, 'range')
df = df.join(colids, 'barcode')
df = df.sort('rowid', 'colid')

#####################
## Export the TSSs ##
#####################

# Create output directory.
output_dir = os.path.splitext(os.path.basename(args.bam))[0]
os.mkdir(output_dir)

# Save the barcodes.
print(f"{datetime.now()} - Saving the barcodes")
colids.drop('colid').coalesce(1).write.text('barcodes')

barcodes_file = [x for x in os.listdir('barcodes') if x.endswith('.txt')][0]
os.rename(f"barcodes/{barcodes_file}", f"{output_dir}/barcodes.tsv")

# Save the features.
print(f"{datetime.now()} - Saving the features")
rowids.drop('rowid').coalesce(1).write.text('features')

features_file = [x for x in os.listdir('features') if x.endswith('.txt')][0]
os.rename(f"features/{features_file}", f"{output_dir}/features.tsv")

# Save the matrix.
print(f"{datetime.now()} - Saving the matrix")
df.select('rowid', 'colid', 'count').coalesce(1).write.option('sep', ' ').csv('matrix')

matrix_file = [x for x in os.listdir('matrix') if x.endswith('.csv')][0]
with open(f"matrix/{matrix_file}", "r+") as mf:
    content = mf.read()
    with open(f"{output_dir}/matrix.mtx", "w") as mf2:
        mf2.write("%%MatrixMarket matrix coordinate integer general\n%\n")
        mf2.write(f"{rowids.count()} {colids.count()} {df.count()}\n")
        mf2.write(content)
