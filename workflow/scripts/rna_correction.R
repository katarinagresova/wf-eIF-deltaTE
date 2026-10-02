# The RNA-seq library correction (config rna_correction:, workflow/rules/rna_correction.smk), on one experiment's RNA-seq
# libraries. Per condition, each library's log2 CPM deviation from the mean of its replicates is fitted as an additive
# smooth of log2 transcript length and GC (median regression on natural splines, on the transcripts with a mean count
# >= min_mean over the condition; predictors clamped to the fitted range, since natural splines extrapolate linearly).
# The fit is the factor a correction would divide out, rescaled so the library keeps its total. The library whose
# factor varies most (SD of its log2 over the fitted transcripts) is corrected if that SD is above max_sd; the others
# are then refitted against the corrected counts, until none is above it. A corrected library's NumReads are divided
# by its factor and its TPM recomputed from them (NumReads / EffectiveLength, scaled to its old TPM total over the same
# rows); every other file is copied unchanged.
# Rscript rna_correction.R <features.tsv> <max_sd> <min_mean> <df> <scan.tsv> <factors.tsv> <n> <in 1..n> <out 1..n>
#                           <condition 1..n>
#   features.tsv  Name, length, gc per transcript
#   out i         <library>_quant.sf: the library names come from these
#   condition i   library i's condition; its replicates are the other libraries of that condition (at least 3 each)
#   scan.tsv      one row per fit: condition, round, library, transcripts fitted, SD of the log2 factor, Spearman of
#                 the deviation with log2 length and with GC before and after the fit, and whether this round
#                 corrected that library
#   factors.tsv   Name + one column per corrected library: the factor its NumReads were divided by
suppressPackageStartupMessages({ library(quantreg); library(splines) })
args <- commandArgs(TRUE)
feat <- read.delim(args[1], row.names = 1)
max_sd <- as.numeric(args[2]); min_mean <- as.numeric(args[3]); df <- as.integer(args[4])
scan_out <- args[5]; factors_out <- args[6]; n <- as.integer(args[7])
stopifnot(length(args) == 7 + 3 * n)
ins <- args[8:(7 + n)]; outs <- args[(8 + n):(7 + 2 * n)]; cond <- args[(8 + 2 * n):(7 + 3 * n)]
libs <- sub("_quant\\.sf$", "", basename(outs))
stopifnot(!anyDuplicated(libs))
q <- setNames(lapply(ins, read.delim, check.names = FALSE), libs)
tx <- q[[1]]$Name
stopifnot(all(sapply(q, function(x) identical(x$Name, tx))), all(tx %in% rownames(feat)))
cnt <- sapply(q, `[[`, "NumReads")
rownames(cnt) <- tx
L <- log2(feat[tx, "length"]); G <- feat[tx, "gc"]
stopifnot(all(table(cond) >= 3))
clamp <- function(x, r) pmin(pmax(x, min(r)), max(r))
rho <- function(a, b) cor(a, b, method = "spearman")

fit_factor <- function(target, refs, cur) {
  m <- cur[, c(target, refs)]
  lcpm <- log2(sweep(m + 0.5, 2, colSums(m), "/") * 1e6)
  rows <- rowMeans(m) >= min_mean
  d <- data.frame(dev = lcpm[rows, target] - rowMeans(lcpm[rows, refs, drop = FALSE]), L = L[rows], G = G[rows])
  fit <- rq(dev ~ ns(L, df = df) + ns(G, df = df), tau = 0.5, data = d)
  f <- 2^predict(fit, newdata = data.frame(L = clamp(L, d$L), G = clamp(G, d$G)))
  f <- f * sum(cur[, target] / f) / sum(cur[, target])  # the library keeps its total
  r <- d$dev - fitted(fit)
  list(f = f, row = data.frame(library = target, transcripts = nrow(d), sd_log2_factor = sd(log2(f[rows])),
                               spearman_length_before = rho(d$dev, d$L), spearman_length_after = rho(r, d$L),
                               spearman_gc_before = rho(d$dev, d$G), spearman_gc_after = rho(r, d$G)))
}

cur <- cnt
scan <- list()
factors <- list()
for (cn in unique(cond)) {
  members <- libs[cond == cn]
  done <- character(0)
  round <- 0
  while (length(done) < length(members)) {
    round <- round + 1
    cand <- setdiff(members, done)
    fits <- setNames(lapply(cand, function(l) fit_factor(l, setdiff(members, l), cur)), cand)
    rows <- do.call(rbind, lapply(fits, `[[`, "row"))
    worst <- rows$library[which.max(rows$sd_log2_factor)]
    fix <- max(rows$sd_log2_factor) > max_sd
    scan[[length(scan) + 1]] <- data.frame(condition = cn, round = round, rows, corrected = fix & rows$library == worst)
    if (!fix) break
    cur[, worst] <- cur[, worst] / fits[[worst]]$f
    factors[[worst]] <- fits[[worst]]$f
    done <- c(done, worst)
  }
}
scan <- do.call(rbind, scan)
print(scan, digits = 3, row.names = FALSE)

for (i in seq_along(libs)) {
  l <- libs[i]
  if (is.null(factors[[l]])) {
    stopifnot(file.copy(ins[i], outs[i], overwrite = TRUE))
    next
  }
  x <- q[[l]]
  stopifnot(all(x$EffectiveLength > 0))
  x$NumReads <- cur[, l]
  rate <- x$NumReads / x$EffectiveLength
  x$TPM <- rate / sum(rate) * sum(q[[l]]$TPM)
  write.table(x, outs[i], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("%s corrected: log2 factor q01 %.2f, q10 %.2f, median %.2f, q90 %.2f, q99 %.2f\n", l,
              quantile(log2(factors[[l]]), .01), quantile(log2(factors[[l]]), .1), median(log2(factors[[l]])),
              quantile(log2(factors[[l]]), .9), quantile(log2(factors[[l]]), .99)))
}
write.table(scan, scan_out, sep = "\t", quote = FALSE, row.names = FALSE)
fac <- data.frame(Name = tx)
for (l in names(factors)) fac[[l]] <- factors[[l]]  # none: Name only
write.table(fac, factors_out, sep = "\t", quote = FALSE, row.names = FALSE)
