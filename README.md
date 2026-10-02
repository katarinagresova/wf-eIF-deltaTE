# wf-eIF-deltaTE

Translational-efficiency (deltaTE) analysis of human Ribo-seq and RNA-seq: per experiment, the
change of every transcript between two conditions in Ribo-seq, in RNA-seq and in translational
efficiency. The Ribo-seq libraries carry a spike-in (yeast in the config and output names), from
which their size factors can be taken.

fastq -> [wf-riboseq-align](https://github.com/katarinagresova/wf-riboseq-align) (trimming,
contaminant filter, STAR, salmon), imported as a Snakemake module -> RNA-seq library correction
-> [ribokit](https://github.com/katarinagresova/ribokit) (Ribo-seq reads per CDS) -> DESeq2
deltaTE model -> master table.

The experiments, their conditions and their libraries are defined in the sample sheet,
`config/samples.csv`: its `experiment` column splits the libraries into as many experiments as
it names, and config `timepoints` selects which of them run (those at the timepoints listed).
Adding an experiment takes only its rows there.

## Alignment: wf-riboseq-align as a module

`workflow/rules/align.smk` imports wf-riboseq-align once, from GitHub at a pinned commit, with
the `align:` block of `config/config.yaml` as its config and `config/samples.csv` as its sample
sheet; its rules are renamed `align_<rule>`. Each experiment gets its own
`results/align/<experiment>/`, and with `mode: filtered` its Ribo-seq reads are mapped to a
transcriptome filtered by that experiment's own RNA-seq. What does not depend on the experiment
(the contaminant index, the RNA-seq index of the unfiltered transcriptome) is built once, in
`results/align/`.

## RNA-seq library correction

`workflow/rules/rna_correction.smk`, per experiment and condition: each RNA-seq library's
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
  align:
    human_transcriptome_fa: "/data/human_transcriptome.fa"
    human_transcriptome_gtf: "/data/human_transcriptome.gtf"
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
