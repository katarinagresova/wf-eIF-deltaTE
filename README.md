# wf-eIF-deltaTE

Translational-efficiency (deltaTE) analysis of the eIF depletion experiments (HCT116,
auxin-induced depletion of eIF3d, eIF4E, eIF4G1, eIF4G2 and eIF4G3, at 4h and 8h; Ribo-seq
with a yeast spike-in, and RNA-seq).

This first part aligns them: fastq ->
[wf-riboseq-align](https://github.com/katarinagresova/wf-riboseq-align) (trimming, contaminant
filter, STAR, salmon), imported as a Snakemake module.

## Alignment: wf-riboseq-align as a module

`workflow/rules/align.smk` imports wf-riboseq-align once, from GitHub at a pinned commit, with
the `align:` block of `config/config.yaml` as its config and `config/samples.csv` as its sample
sheet; its rules are renamed `align_<rule>`. The sample sheet's `experiment` column
(`<target>_<timepoint>`) splits the libraries into the 10 experiments: each gets its own
`results/align/<experiment>/`, and with `mode: filtered` its Ribo-seq reads are mapped to a
transcriptome filtered by that experiment's own RNA-seq. What does not depend on the experiment
(the contaminant index, the RNA-seq index of the unfiltered transcriptome) is built once, in
`results/align/`.

Adding an experiment takes only its rows in `config/samples.csv`.

## Setup

- The paths in `config/config.yaml` and `config/samples.csv` are `/path/to/` placeholders. Set
  them there, or put the real ones in an untracked `config/local/config.yaml`, which
  `snakemake.sh` and `slurm_job.sh` pass with `--configfile` when it exists. It needs only the
  keys it changes:
  ```yaml
  samples: "config/local/samples.csv"
  align:
    contaminants_fa: "/data/contaminants.fa"
    human_transcriptome_fa: "/data/human_transcriptome.fa"
  ```
- The contaminant set (`align: contaminants_fa`) is not built by the run. wf-riboseq-align ships
  one, `resources/contaminants_built.fa`; to rebuild it from public sources, run its
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
