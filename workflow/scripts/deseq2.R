# DESeq2 on one experiment's Ribo-seq and RNA-seq libraries (workflow/rules/deltate.smk), from deltate_inputs.py's
# tables. Fold changes are treatment over reference. Five models:
#   deltaTE_<norm>       the change in translational efficiency (Chothani et al. 2019): both assays together, design
#                        ~ condition + assay + condition:assay, the coefficient condition<treatment>:assayribo
#   diffribo_<norm>      the Ribo-seq alone, design ~ condition
#   difftotal_autonorm   the RNA-seq alone, design ~ condition
# <norm> is how the Ribo-seq size factors are set: yeastnorm, from the spike-in (each library's spike-in reads over the
# first Ribo-seq library's); autonorm, estimated within the Ribo-seq libraries. The RNA-seq ones are always estimated
# within the RNA-seq libraries.
# Estimated within an assay = as estimateSizeFactors() does for tximport input: median of ratios on counts / length,
# the lengths row-centred, and the model gets size factor x row-centred length as normalization factors. The deltaTE
# models get them as such (lengths row-centred over all libraries); diffribo_yeastnorm gets plain size factors (the
# Ribo-seq length is the CDS's, the same in every Ribo-seq library).
# The Ribo-seq is counted on CDSs (a transcript without one is NA in its columns): the deltaTE and Ribo-seq models are
# fitted on the transcripts with a CDS, the RNA-seq model and size factors on every transcript.
# results() tunes its independent filtering to alpha = 0.05, the padj cutoff of a hit (DESeq2's default: 0.1).
# Rscript deseq2.R <dir> <reference> <treatment>
# Reads <dir>/{sampleTable,counts,length,spike_in_reads}.tsv; writes <dir>/deseq_res_<model>.tsv and
# <dir>/plots/MAplot_diff_{deltaTE,ribo,total}_<norm>.pdf.
suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
})

args <- commandArgs(TRUE)
stopifnot(length(args) == 3)
dir <- args[1]; conditions <- args[2:3]

read_matrix <- function(file) {
  x <- read_tsv(file.path(dir, file), col_types = cols(Name = "c", .default = "d"))
  m <- as.matrix(x[-1]); rownames(m) <- x$Name
  m
}
samples <- as.data.frame(read_tsv(file.path(dir, "sampleTable.tsv"), col_types = "ccc"))
rownames(samples) <- samples$sampleName
samples$assay <- factor(samples$assay, levels = c("total", "ribo"))
samples$condition <- factor(samples$condition, levels = conditions)
counts <- round(read_matrix("counts.tsv")); mode(counts) <- "integer"
lengths <- read_matrix("length.tsv")
spike_in <- read_tsv(file.path(dir, "spike_in_reads.tsv"), col_types = "cd")
ribo <- samples$sampleName[samples$assay == "ribo"]
total <- samples$sampleName[samples$assay == "total"]
stopifnot(!anyNA(samples), identical(colnames(counts), samples$sampleName),
          identical(dimnames(lengths), dimnames(counts)), identical(is.na(lengths), is.na(counts)),
          identical(spike_in$sampleName, ribo))
cds <- complete.cases(counts)

# a DESeqDataSet as DESeqDataSetFromTximport() builds it
dataset <- function(rows, libs, design) {
  dds <- DESeqDataSetFromMatrix(counts[rows, libs], samples[libs, ], design)
  assay(dds, "avgTxLength") <- lengths[rows, libs]
  dds
}
centre <- function(len) len / exp(rowMeans(log(len)))
alpha <- 0.05

save <- function(res, model, plot) {
  message(model, ": ", sum(res$padj < alpha, na.rm = TRUE), " of ", nrow(res), " at padj < ", alpha)
  write_tsv(data.frame(Name = rownames(res), as.data.frame(res)), file.path(dir, paste0("deseq_res_", model, ".tsv")))
  pdf(file.path(dir, "plots", paste0("MAplot_diff_", plot, ".pdf")), height = 6, width = 8)
  plotMA(res, ylim = c(-8, 8), colSig = "red")
  invisible(dev.off())
}
dir.create(file.path(dir, "plots"), showWarnings = FALSE)

sf_total <- estimateSizeFactorsForMatrix(counts[, total] / centre(lengths[, total]))
sf_ribo <- list(yeastnorm = setNames(spike_in$reads / spike_in$reads[1], spike_in$sampleName),
                autonorm = estimateSizeFactorsForMatrix(counts[cds, ribo] / centre(lengths[cds, ribo])))
print(list(total = sf_total, ribo = sf_ribo))

te <- dataset(cds, samples$sampleName, ~ condition + assay + condition:assay)
te_coef <- make.names(paste0("condition", conditions[2], ":assayribo"))
for (norm in names(sf_ribo)) {
  sf <- c(sf_ribo[[norm]], sf_total)
  stopifnot(identical(names(sf), colnames(te)))
  normalizationFactors(te) <- sweep(centre(lengths[cds, ]), 2, sf, "*")
  save(results(DESeq(te), name = te_coef, alpha = alpha), paste0("deltaTE_", norm), paste0("deltaTE_", norm))
}

save(results(DESeq(dataset(rownames(counts), total, ~ condition)), alpha = alpha),
     "difftotal_autonorm", "total_autonorm")
save(results(DESeq(dataset(cds, ribo, ~ condition)), alpha = alpha), "diffribo_autonorm", "ribo_autonorm")
ribo_yeast <- dataset(cds, ribo, ~ condition)
sizeFactors(ribo_yeast) <- sf_ribo$yeastnorm
save(results(DESeq(ribo_yeast), alpha = alpha), "diffribo_yeastnorm", "ribo_yeastnorm")
