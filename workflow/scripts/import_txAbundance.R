# One experiment's Ribo-seq and RNA-seq quantifications as one tximport object (workflow/rules/deltate.smk), Ribo-seq
# libraries first, on the transcripts of the first RNA-seq library (a transcript a Ribo-seq quant lacks gets 0 reads
# and length 1), and the sample table DESeq2 models them with, split from the sample ids:
# <assay>_<id>_<subunit>_<auxin>_<harringtonine>_<time>_<replicate> (deltate.smk checks every id).
# Rscript import_txAbundance.R <dir> <n_ribo> <n_rna> <ribo sample 1..n_ribo> <rna sample 1..n_rna>
#                              <ribo quant 1..n_ribo> <rna quant 1..n_rna>
#   ribo sample i  sample id, ending in "_R1" (dropped)
#   ribo quant i   reads per human CDS (columns Name Length EffectiveLength ritpm NumReads)
#   rna quant i    salmon's quant.sf
# Writes into <dir>: sampleTable.{rds,tsv}, txi.rds, and its matrices as txi_{counts,abundance,length}.tsv.
library(magrittr)
library(tximport)
library(tidyverse)

args <- commandArgs(TRUE)
out_dir <- args[1]; n_ribo <- as.integer(args[2]); n_rna <- as.integer(args[3])
stopifnot(length(args) == 3 + 2 * (n_ribo + n_rna))
ribo_samples <- args[4:(3 + n_ribo)]
rna_samples <- args[(4 + n_ribo):(3 + n_ribo + n_rna)]
ribo_files <- args[(4 + n_ribo + n_rna):(3 + 2 * n_ribo + n_rna)]
rna_files <- args[(4 + 2 * n_ribo + n_rna):length(args)]


# Read salmon quantification output .tsv file
read_quantfile_rna <- function(filepath){
  read_tsv(filepath) %>%
    filter(Name %in% my_transcripts) %>%
    arrange(Name) %>%
    mutate(TPM      = replace_na(TPM, 0),
           NumReads = replace_na(NumReads, 0))
}

# Import ribo and RNA separately, using the full RNA transcript universe obtained from salmon output files
read_quantfile_ribo <- function(filepath){
  read_tsv(filepath) %>%
    # Right-join to full RNA universe so missing transcripts get NA → 0
    right_join(tibble(Name = my_transcripts), by = "Name") %>%
    arrange(Name) %>%
    rename_at(vars(contains("ritpm")), list(~ sub("ritpm", "TPM", .))) %>%
    mutate(TPM      = replace_na(TPM, 0),
           NumReads = replace_na(NumReads, 0),
           Length   = replace_na(Length, 1),
           EffectiveLength = replace_na(EffectiveLength, 1))   # DESeq2 needs a non-zero length in the EffectiveLength column
}


# construct sample metadata table
sampleCols <- c("assay", "sample_id", "subunit", "auxin", "harringtonine", "time", "replicate")
ribo_samples <- str_sub(ribo_samples, 1, -4)
samples <- c(ribo_samples, rna_samples)
sampleTable <- data.frame(samples) %>%
  rename(sampleName = samples) %>%
  separate(col = sampleName, into = sampleCols, sep = "_", remove = FALSE)
print(sampleTable)

#define transcript universe based on first sample (they all share the same transcriptome)
my_transcripts <- read_tsv(rna_files[1]) %>% .$Name


# name the ribo and salmon quantifications by sample
names(ribo_files) <- ribo_samples
names(rna_files) <- rna_samples
quantfiles <- c(ribo_files, rna_files)
print(quantfiles)

# create txi objects
txi_ribo <- tximport(ribo_files,  type = "salmon", txOut = TRUE, importer = read_quantfile_ribo)
txi_rna  <- tximport(rna_files,   type = "salmon", txOut = TRUE, importer = read_quantfile_rna)
# Recombine into one txi — columns are already ordered by Name because of arrange()
txi <- list(
  counts    = cbind(txi_ribo$counts,    txi_rna$counts),
  abundance = cbind(txi_ribo$abundance, txi_rna$abundance),
  length    = cbind(txi_ribo$length,    txi_rna$length),
  countsFromAbundance = txi_rna$countsFromAbundance
)
class(txi) <- "list"
print(colnames(txi$counts))
print(dim(txi$abundance))


# Save sampleTable
saveRDS(sampleTable, file.path(out_dir, "sampleTable.rds"))
write_tsv(sampleTable, file.path(out_dir, "sampleTable.tsv"))

# Save txi matrices — each is a transcript=rows x sample=cols matrix so need to convert to a data frame
saveRDS(txi, file.path(out_dir, "txi.rds"))
as.data.frame(txi$counts)    %>% rownames_to_column("Name") %>% write_tsv(file.path(out_dir, "txi_counts.tsv"))
as.data.frame(txi$abundance) %>% rownames_to_column("Name") %>% write_tsv(file.path(out_dir, "txi_abundance.tsv"))
as.data.frame(txi$length)    %>% rownames_to_column("Name") %>% write_tsv(file.path(out_dir, "txi_length.tsv"))
