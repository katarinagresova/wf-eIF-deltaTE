# The change in translational efficiency (workflow/rules/deltate.smk): DESeq2 on the Ribo-seq and RNA-seq libraries
# together, design ~ auxin + assay + auxin:assay, the interaction auxinplusAux.assayribo being the change in TE. Twice:
# the Ribo-seq size factors from the spike-in (yeastnorm, calc_sizefactors.R's), and estimated within the Ribo-seq
# libraries (autonorm); the RNA-seq size factors are estimated within the RNA-seq libraries both times.
# Rscript run_deltaTE.R <dir>
# Reads <dir>/sampleTable.rds, txi.rds (import_txAbundance.R) and sizeFactors_yeast_ribo.rds; writes into <dir>:
# deseq_res_deltaTE_{yeastnorm,autonorm}.{tsv,rds}, sizeFactors_human_{total,ribo_auto}.rds, and
# plots/MAplot_diff_deltaTE_{yeastnorm,autonorm}.pdf.

# the workflow employed here is based on Chotani et. al 2019 Current Protocols in Molecular Biology

library(magrittr)
library(DESeq2)
library(tidyverse)

out_dir <- commandArgs(TRUE)[1]

#read sample metadata table created with import_txAbundance.R
sampleTable <- readRDS(file.path(out_dir, "sampleTable.rds"))
print(sampleTable)
#read txi object created with import_txAbundance.R
txi <- readRDS(file.path(out_dir, "txi.rds"))
#print(txi)

#build DESeq dataset
dds <- DESeqDataSetFromTximport(txi = txi, colData = sampleTable, design = ~ auxin+assay+auxin:assay)
#change reference level to RNA-seq (total rna)
dds$assay <- relevel(dds$assay, ref = "total")
# The size factors below are estimated on counts / avgTxLength, as DESeq2 does for tximport input, so they go in
# together with avgTxLength, as normalization factors (size factor x avgTxLength, row-centred), which is what DESeq2's
# own estimateSizeFactors() builds; set with sizeFactors<-, DESeq() would apply them to the counts without the lengths.
# The Ribo-seq lengths are each CDS's length (or 1), the same in every Ribo-seq library, so there it changes nothing.
normMatrix_all <- assays(dds)[["avgTxLength"]]
normMatrix_all <- normMatrix_all / exp(rowMeans(log(normMatrix_all)))
with_lengths <- function(sf) {
  stopifnot(identical(names(sf), colnames(normMatrix_all)))
  sweep(normMatrix_all, 2, sf, "*")
}


# Make txi and sample table with only RNA samples
sampleTable_rnaseq <- sampleTable %>% filter(assay == "total")
txi_rnaseq <- txi
keep <- grepl("total", colnames(txi$counts))
txi_rnaseq$counts <- txi$counts[, keep, drop = FALSE]
txi_rnaseq$abundance <- txi$abundance[, keep, drop = FALSE]
txi_rnaseq$length <- txi$length[, keep, drop = FALSE]

# Build DESeq data set for just RNA
dds_rnaseq <- DESeqDataSetFromTximport(txi = txi_rnaseq, colData = sampleTable_rnaseq, design = ~ auxin)

# Run DESeq2 on RNA only to extract sizefactors and save
dds_rnaseq <- DESeq(dds_rnaseq)
normMatrix <- assays(dds_rnaseq)[["avgTxLength"]]
normMatrix <- normMatrix / exp(rowMeans(log(normMatrix)))
my_sizeFactors_rna <- estimateSizeFactorsForMatrix(counts(dds_rnaseq) / normMatrix)
saveRDS(my_sizeFactors_rna, file.path(out_dir, "sizeFactors_human_total.rds"))

# Retrieve RPF size factors created with calc_sizefactors.R
ribo_samples <- sampleTable %>% filter(assay == "ribo") %>% .$sampleName
my_sizeFactors_ribo <- readRDS(file.path(out_dir, "sizeFactors_yeast_ribo.rds"))
my_sizeFactors_ribo <- my_sizeFactors_ribo[ribo_samples]

# Manually inject size factors to DESeqDataSet
dds_manual <- dds
normalizationFactors(dds_manual) <- with_lengths(c(my_sizeFactors_ribo, my_sizeFactors_rna))

# Run DESeq2 and manually extract statistics for interaction term, which represents deltaTE
dds_manual <- DESeq(dds_manual)
res_manual <- results(dds_manual, name = "auxinplusAux.assayribo") #this step is not possible in python :()
print("Dimensions of results table with manual sizefactor injection: transcripts(=rows) x DEseq2 statistics(=columns)")
print(dim(res_manual))

# save results to .rds and .tsv
saveRDS(res_manual, file.path(out_dir, "deseq_res_deltaTE_yeastnorm.rds"))
res_manual_df <- as.data.frame(res_manual) %>% rownames_to_column(var = "Name")
write_tsv(res_manual_df, file = file.path(out_dir, "deseq_res_deltaTE_yeastnorm.tsv"))

#plotting MA plot
plotMAplot <- function(INPUT, ylim, filename){
  plotfile <- file.path(out_dir, "plots", filename) %T>% pdf(h=6, w=8)

  rwplot <- plotMA(INPUT, ylim = ylim, colSig = "red")

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}

res_manual %>% plotMAplot(., ylim = c(-8, 8), filename = "MAplot_diff_deltaTE_yeastnorm.pdf")


# DESeq2 with automatic size factors (autonorm): within-ribo length-corrected size factors ---------
# Replaces the yeast spike-in for RPF. RNA size factors (my_sizeFactors_rna) are reused from above.

# Make txi and sample table with only RPF samples (mirror of the RNA-only block above)
txi_ribo <- txi
keep_ribo <- grepl("ribo", colnames(txi$counts))
txi_ribo$counts <- txi$counts[, keep_ribo, drop = FALSE]
txi_ribo$abundance <- txi$abundance[, keep_ribo, drop = FALSE]
txi_ribo$length <- txi$length[, keep_ribo, drop = FALSE]
sampleTable_ribo <- sampleTable %>% filter(assay == "ribo")

# Build RPF-only DESeq dataset and estimate length-corrected size factors within RPF only
dds_ribo <- DESeqDataSetFromTximport(txi = txi_ribo, colData = sampleTable_ribo, design = ~ auxin)
dds_ribo <- DESeq(dds_ribo)
normMatrix_ribo <- assays(dds_ribo)[["avgTxLength"]]
normMatrix_ribo <- normMatrix_ribo / exp(rowMeans(log(normMatrix_ribo)))
my_sizeFactors_ribo_auto <- estimateSizeFactorsForMatrix(counts(dds_ribo) / normMatrix_ribo)
my_sizeFactors_ribo_auto <- my_sizeFactors_ribo_auto[ribo_samples]
saveRDS(my_sizeFactors_ribo_auto, file.path(out_dir, "sizeFactors_human_ribo_auto.rds"))

# Inject automatic size factors (RPF then RNA, matching column order) and run DESeq2
dds_auto <- dds
normalizationFactors(dds_auto) <- with_lengths(c(my_sizeFactors_ribo_auto, my_sizeFactors_rna))
dds_auto <- DESeq(dds_auto)
res_auto <- results(dds_auto, name = "auxinplusAux.assayribo")
print("Dimensions of results table with automatic sizefactors: transcripts(=rows) x DEseq2 statistics(=columns)")
print(dim(res_auto))

# save results to .rds and .tsv
saveRDS(res_auto, file.path(out_dir, "deseq_res_deltaTE_autonorm.rds"))
res_auto_df <- as.data.frame(res_auto) %>% rownames_to_column(var = "Name")
write_tsv(res_auto_df, file = file.path(out_dir, "deseq_res_deltaTE_autonorm.tsv"))

res_auto %>% plotMAplot(., ylim = c(-8, 8), filename = "MAplot_diff_deltaTE_autonorm.pdf")
