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
# First argument is output report name, second is CRISPResso2 dir, third is 3ko vs perv report
sample_bgh <- args[1]
report_name <- strsplit(basename(sample_bgh), "_")[[1]][1]
report_type <- "payload_aberrant_insertions"

# Ensure CRISPResso2 output exists
if (!file.exists(sample_bgh)) {
  print(str_glue("Couldn't find {sample_bgh} "))
  quit(status = 1)
}

report_name <- paste0(report_name, "_", report_type, "_report")
rmarkdown::render("coverage_bedgraph.Rmd",
                  params = list(sample_bedgraph = sample_bgh),
                  output_file = report_name)
