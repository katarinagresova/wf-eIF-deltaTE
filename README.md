# wf-eIF-deltaTE

Translational-efficiency (deltaTE) analysis of human Ribo-seq and RNA-seq: per experiment, the
change of every transcript between two conditions in Ribo-seq, in RNA-seq and in translational
efficiency. The Ribo-seq libraries carry a spike-in (yeast in the config and output names), from
which their size factors can be taken.

Ribo-seq fastq -> [wf-riboseq-align](https://github.com/katarinagresova/wf-riboseq-align)
(trimming, contaminant filter, bowtie), imported as a Snakemake module ->
[ribokit](https://github.com/katarinagresova/ribokit) (Ribo-seq reads per CDS); RNA-seq fastq ->
salmon -> expression filter -> RNA-seq library correction; both -> DESeq2 deltaTE model ->
master table. Both assays are counted on the same transcriptome, collapsed by coding sequence.

The experiments, their conditions and their libraries are defined in the sample sheet,
`config/samples.csv`: its `experiment` column splits the libraries into as many experiments as
it names, and config `timepoints` selects which of them run (those at the timepoints listed).
Adding an experiment takes only its rows there.

## Reference transcriptome

`workflow/rules/collapse_transcriptome.smk` collapses the human transcriptome of config
`collapse_transcriptome` (fasta, GTF, and an external short-read TPM table for the tie-break).
First, the genes that GENCODE (`gencode_gtf`) calls readthroughs are dropped (RPS10-NUDT3): their
CDS spans two genes' CDSs. Then transcripts with an identical CDS keep one, across gene ids too, as a
Ribo-seq read cannot tell them apart; a CDS at least 99% covered by another of the same gene merges
into it; and each gene keeps one transcript. Last, a CDS of which fewer than `min_unique_fraction`
(0.2) of the `kmer`-mers (28) occur in no other CDS merges into the CDS it shares most with
(paralogs: EIF3C / EIF3CL, SMN1 / SMN2). A merged group keeps a protein-coding gene, then the higher
TPM. Readthroughs and paralogs share most of their footprints with another gene, so their counts
jumped between replicates. `results/resources/` holds the collapsed fasta and GTF, and
`collapse_report.tsv` (every group, and why it was resolved). Everything downstream uses it.

## Alignment: wf-riboseq-align as a module

`workflow/rules/align.smk` imports wf-riboseq-align once, from GitHub at a pinned commit, with
the `align:` block of `config/config.yaml` as its config and `config/samples.csv` as its sample
sheet; its rules are renamed `align_<rule>`. It aligns the Ribo-seq only (the `rna` rows are
skipped), to the whole collapsed transcriptome plus the spike-in, and splits the alignments into a
human and a spike-in BAM file. Each experiment gets its own `results/align/<experiment>/` (split
BAMs, `qc/summary.tsv`, FastQC, MultiQC). What does not depend on the experiment (the Ribo-seq
reference and the combined bowtie index of the transcriptome, spike-in and contaminants) is built
once, in `results/align/`.

## RNA-seq: salmon and the expression filter

`workflow/rules/rnaseq.smk`, per experiment (the rules were wf-riboseq-align's until it became
Ribo-seq only, behaviour unchanged):
- salmon 1.10.2 (the version of the original pipeline's container) quantifies each paired-end
  library on its fastq files, by selective alignment (`-l A --seqBias --gcBias`), against an index
  of the collapsed transcriptome with every sequence of the genome (config
  `references: human_genome_fa`) as a decoy: a fragment that aligns better to the genome than to
  any transcript is counted for none. `--keepDuplicates`, so that a transcript whose sequence
  duplicates another's stays in `quant.sf`. salmon samples fragments for its bias models from a
  random seed it does not expose, so a rerun can differ slightly (at most a few transcripts change
  sides of the expression filter's threshold); `-p 1` does not make it deterministic either.
- The expression filter (config `autofilter`): a transcript is kept if its TPM is at least
  `min_tpm` in at least `min_samples` of the experiment's RNA-seq libraries, else blacklisted.
  The blacklist applies to what is counted, not to the Ribo-seq alignment: ribokit counts on the
  collapsed GTF minus the blacklist, and the RNA-seq library correction, then DESeq2, read
  `quant.rnaseq_filtered.sf`: salmon's `quant.sf` minus the blacklisted rows, sorted by name, the
  values as salmon wrote them (TPM is not renormalised).

The index is built once, in `results/rnaseq/salmon_index/`; `results/rnaseq/<experiment>/` holds
`salmon/<sample>/` (salmon's output and `quant.rnaseq_filtered.sf`),
`rnaseq_filter_blacklist_txid.txt`, `human_transcriptome.rnaseq_filtered.gtf`, FastQC of the raw
fastq files (`fastqc/`) and a MultiQC report of salmon and FastQC (`multiqc/`).

## RNA-seq library correction

`workflow/rules/rna_correction.smk`, per experiment and condition, on the filtered quants: each
RNA-seq library's
deviation from its replicates is fitted as a smooth function of transcript length and GC. While
the most variable fit has an SD above `max_sd`, that library is divided by its fit and the
others are refitted (config `rna_correction`). `results/rna_correction/<experiment>/` holds
every library's quant (corrected or copied as is), `scan.tsv` (every fit) and `factors.tsv`.

## Ribo-seq quantification: ribokit

`workflow/rules/ribokit.smk`, per Ribo-seq library: ribokit estimates the P-site offsets from
the reads on the human transcripts and counts the reads per CDS, human and spike-in (with the
human offsets). Only reads of config `read_lengths` count, the same window in every library.
Reads compatible with several CDSs are split by an EM that starts deterministically and fails if
it does not converge. `results/ribokit/<experiment>/` holds
`<sample>.{human,yeast}.{quant,offsets,ties,stats,psites}.tsv` (`psites.tsv`: one row per
alignment assigned to a CDS, its P-site position). ribokit is installed from GitHub at the
commit pinned in `workflow/envs/ribokit.yaml`.

## deltaTE: DESeq2

`workflow/rules/deltate.smk`, per experiment: the deltaTE model (Chotani et al. 2019), DESeq2 on
the Ribo-seq and RNA-seq libraries together with design `~ condition + assay + condition:assay`;
the change in translational efficiency is the interaction. Fold changes are treatment over
reference (config `conditions`). The Ribo-seq size factors come from the spike-in (`yeastnorm`)
or are estimated within the libraries (`autonorm`), and go into the model together with the
transcript lengths, as normalization factors. Each assay is also fitted alone.
`results/deltaTE/<experiment>/` holds (as `.tsv` and `.rds`, MA plots in `plots/`):

| table | what |
|---|---|
| `deseq_res_deltaTE_{yeastnorm,autonorm}` | change in translational efficiency (log2FC_TE) |
| `deseq_res_diffribo_{yeastnorm,autonorm}` | Ribo-seq alone (log2FC_RPF) |
| `deseq_res_difftotal_autonorm` | RNA-seq alone (log2FC_RNA) |

## Master table

`workflow/rules/master.smk` -> `results/master/eif-master.csv`: all experiments run, per
transcript: log2FC and padj of TE (`deltaTE_yeastnorm`), RPF (`diffribo_yeastnorm`) and RNA
(`difftotal_autonorm`), the TPM of every library (ribokit's ritpm, salmon's TPM), log2TE per
replicate, and the transcript annotation from the GTF. The experiments are stacked by target
(the targets in `MASTER_TARGET_ORDER` in master.smk first, any other by name) and timepoint.

`results/master/eif-master.filtered.csv`: the same table restricted, per experiment, to the
transcripts with Ribo-seq ritpm >= 5 and RNA-seq TPM >= 1 in every library of the reference
condition (config `master_filter`). The Ribo-seq floor is the paper's row filter; the RNA-seq
floor drops transcripts whose TE is meaningless because the RNA-seq barely sees them (mostly
histone mRNAs). Filtered after DESeq2, so padj is the full table's.

## Sample sheet

`config/samples.csv` has one row per library: `sample_id`, `target`, `timepoint`, `condition`,
`rep`, `read_type` (`ribo` or `rna`), `fastq_1`, `fastq_2` (empty for single-end) and
`experiment`.
- Experiment names are free (e.g. `<target>_<timepoint>`); nothing parses them, nor the sample
  ids. Assay, condition, target and timepoint come from their columns.
- Every experiment run needs both conditions of config `conditions`, and no other, in its
  Ribo-seq and in its RNA-seq libraries (checked when the workflow starts).
- Row order matters: the spike-in size factors scale every Ribo-seq library to the experiment's
  first Ribo-seq row.

## Setup

- The paths in `config/config.yaml` and `config/samples.csv` are `/path/to/` placeholders. Set
  them there, or put the real ones in an untracked `config/local/config.yaml`, which
  `snakemake.sh` and `slurm_job.sh` pass with `--configfile` when it exists. It needs only the
  keys it changes:
  ```yaml
  samples: "config/local/samples.csv"
  collapse_transcriptome:
    raw_fa: "/data/human_transcriptome.fa"
    raw_gtf: "/data/human_transcriptome.gtf"
    tpm_table: "/data/transcript_tpm_short_read.csv"
    gencode_gtf: "/data/gencode.v47.annotation.gtf"
  align:
    spike_in_transcriptome_fa: "/data/yeast_transcriptome.fa"
  references:
    human_genome_fa: "/data/human_genome.fa"
    yeast_genome_fa: "/data/yeast_genome.fa"
    yeast_gtf: "/data/yeast.gtf"
  ```
- The contaminant set is not built by the run: it is the one wf-riboseq-align ships,
  `resources/contaminants_built.fa`, fetched from GitHub at the pinned commit.
  `align: contaminants_fa: <fasta>` replaces it; to rebuild it from public sources, run its
  `build_contaminants` rule (needs internet access), e.g.
  `./snakemake.sh results/align/reference/contaminants_built.fa`.
- `align_snakefile: <checkout>/workflow/Snakefile` imports a local wf-riboseq-align checkout
  instead of the pinned commit.

## Running

Needs a conda env `snake` with snakemake 9 (and snakemake-executor-plugin-slurm for SLURM);
every rule runs in its own conda env.

```
./snakemake.sh -n        # dry run
./snakemake.sh           # local run (workflow/profiles/default)
sbatch slurm_job.sh      # one SLURM job per rule (workflow/profiles/slurm)
```

## License

Apache License 2.0 (`LICENSE`, `NOTICE`).
