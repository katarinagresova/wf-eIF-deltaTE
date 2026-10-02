# Differential expression within one assay (workflow/rules/deltate.smk): DESeq2 on its libraries alone, design
# ~ condition (treatment over reference, the sample table's factor levels), with size factors estimated within them
# (autonorm), and for ribo also from the spike-in (yeastnorm, calc_sizefactors.R's).
# Rscript run_deseq.R <dir> <assay>
#   assay  total (RNA-seq) or ribo (Ribo-seq)
# Reads <dir>/sampleTable.rds, txi.rds (import_txAbundance.R) and, for ribo, sizeFactors_yeast_ribo.rds; writes into
# <dir>: deseq_res_diff<assay>_autonorm.{tsv,rds} (+ _yeastnorm for ribo) and plots/MAplot_diff_<assay>_<norm>.pdf.
library(magrittr)
library(DESeq2)
library(tidyverse)

#parse CLI arguments
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
my_assay <- args[2]
stopifnot(my_assay %in% c("total", "ribo"))


# function to save MA plot
plotMAplot <- function(INPUT, ylim, filename){
  plotfile <- file.path(out_dir, "plots", filename) %T>% pdf(h=6, w=8)

  rwplot <- plotMA(INPUT, ylim = ylim, colSig = "red")

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


# load general sampleTable including both RNA and RPF created by import_txAbundance.R and filter to selected assay type only
sampleTable <- readRDS(file.path(out_dir, "sampleTable.rds"))
sampleTable_my_assay <- sampleTable %>% filter(assay == my_assay)

# load general txi object including both RNA and RPF transcripts created by import_txAbundance.R
# and keep only columns corresponding to samples of selected assay type
txi <- readRDS(file.path(out_dir, "txi.rds"))
txi_my_assay <- txi
keep <- colnames(txi$counts) %in% sampleTable_my_assay$sampleName
txi_my_assay$counts <- txi$counts[, keep, drop = FALSE]
txi_my_assay$abundance <-txi$abundance[, keep, drop = FALSE]
txi_my_assay$length <- txi$length[, keep, drop = FALSE]




# Build DESeq data set
dds <- DESeqDataSetFromTximport(txi = txi_my_assay, colData = sampleTable_my_assay, design = ~ condition)




# DESeq2 with automatic size factor estimation (autonorm) - runs for any assay --------------------
# For a single-assay dataset, DESeq()'s automatic estimation uses the length-corrected
# normalizationFactors from tximport, i.e. within-assay length-corrected normalization (method B).
dds_auto <- DESeq(dds)
res_auto <- results(dds_auto)

# MA plots
res_auto %>% plotMAplot(., ylim = c(-8, 8), filename = paste0("MAplot_diff_", my_assay, "_autonorm.pdf"))

# save results to .rds and .tsv
saveRDS(res_auto, file.path(out_dir, paste0("deseq_res_diff", my_assay, "_autonorm.rds")))
res_auto_df <- as.data.frame(res_auto) %>% rownames_to_column(var = "Name")
write_tsv(res_auto_df, file = file.path(out_dir, paste0("deseq_res_diff", my_assay, "_autonorm.tsv")))


# DESeq2 with manual sizefactors from yeast spike-in for RPF (yeastnorm) - ribo only ---------------
if(my_assay == "ribo"){

# load manual sizefactors from calc_sizefactors.R and filter to RPF only
my_sizeFactors <- readRDS(file.path(out_dir, "sizeFactors_yeast_ribo.rds"))
my_sizeFactors <- my_sizeFactors[sampleTable_my_assay$sampleName]

# Manually set size factors
dds_manual <- dds
sizeFactors(dds_manual) <- my_sizeFactors

# Run DESeq2
dds_manual <- DESeq(dds_manual)
res_manual <- results(dds_manual)

# MA plots
res_manual %>% plotMAplot(., ylim = c(-8, 8), filename = paste0("MAplot_diff_", my_assay, "_yeastnorm.pdf"))

#save results to .rds and .tsv
saveRDS(res_manual, file.path(out_dir, paste0("deseq_res_diff", my_assay, "_yeastnorm.rds")))
res_manual_df <- as.data.frame(res_manual) %>% rownames_to_column(var = "Name")
write_tsv(res_manual_df, file = file.path(out_dir, paste0("deseq_res_diff", my_assay, "_yeastnorm.tsv")))
}
