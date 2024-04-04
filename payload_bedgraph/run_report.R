#!/usr/bin/env Rscript

### run_report.R ###
#
# This script takes in th bedgraph file generated from mapping PL15S reads against Sus scrofa v11.1 to identify aberrant PL15S insertions
# then passes those to generate a dynamic report using R markdown and the template file
# in coverage_bedgraph.Rmd, depending on the
# type of report specified

library(dplyr)
library(rmarkdown)

args <- commandArgs(trailingOnly = TRUE)
# First argument is sample bedgraph file
sample_bgh <- args[1]
min_DP <- args[2]
stat1 <- args[3] # BAM stats after aligning PacBio reads to Payload plasmid
stat2 <- args[4] # 2nd alignment step, BAM stats after aligning payload reads to Ssc11.1, removing homologous regions
report_name <- strsplit(basename(sample_bgh), "_")[[1]][1]
report_type <- "payload_aberrant_insertions"

# Ensure bedgraph exists
if (!file.exists(sample_bgh)) {
  print(str_glue("Couldn't find {sample_bgh} "))
  quit(status = 1)
}

if (!file.exists(stat1)) {
  print(str_glue("Couldn't find {stat1} "))
  quit(status = 1)
}

if (!file.exists(stat2)) {
  print(str_glue("Couldn't find {stat2} "))
  quit(status = 1)
}

report_name <- paste0(report_name, "_", report_type, "_report")
rmarkdown::render("/opt/coverage_bedgraph.Rmd",
                  params = list(sample_bedgraph = sample_bgh, DP_cutoff = min_DP, pl_mapped_stat = stat1, ssc11_mapped_stat = stat2 ),
                  output_file = report_name)
