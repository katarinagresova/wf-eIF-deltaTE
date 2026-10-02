# The Ribo-seq size factors from the spike-in (workflow/rules/deltate.smk): each library's spike-in (yeast) CDS reads
# over the first library's, so the order of the libraries matters.
# Rscript calc_sizefactors.R <dir> <n> <ribo sample 1..n> <yeast quant 1..n>
#   ribo sample i  library i's sample id (the size factors' names)
#   yeast quant i  library i's reads per spike-in CDS (columns Name, NumReads)
# Writes <dir>/sizeFactors_yeast_ribo.rds.
library(tidyverse)

args <- commandArgs(TRUE)
out_dir <- args[1]; n <- as.integer(args[2])
stopifnot(length(args) == 2 + 2 * n)
ribo_samples <- args[3:(2 + n)]
ribo_files <- setNames(args[(3 + n):(2 + 2 * n)], ribo_samples)

# Compile read counts into a list of tables for all yeast RPF samples
ribostan_tables <- list()
for(sample in ribo_samples){
  ribostan_tables[[sample]] <- read.table(ribo_files[[sample]], header = TRUE)
  ribostan_tables[[sample]] <- ribostan_tables[[sample]] %>%
    select(Name, NumReads) # keep only transcript {Name} and number of raw reads {NumReads} columns
}

# Collapse list into a data frame
ribo_readcounts_df <- bind_rows(ribostan_tables, .id = "sample") %>%
  pivot_wider(names_from = sample, values_from = NumReads)
print(dim(ribo_readcounts_df))

# Calculate size factors relative to first library:
# Calculate column-wise (=sample-wise) sums and divide sum of reads for each sample by sum or reads in first sample
sumReads <- colSums(ribo_readcounts_df[,-1], na.rm = TRUE)
sizeFactors <- sumReads/sumReads[1]
print(sizeFactors)
saveRDS(sizeFactors, file.path(out_dir, "sizeFactors_yeast_ribo.rds"))
