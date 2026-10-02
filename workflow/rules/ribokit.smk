# Ribo-seq quantification with ribokit (envs/ribokit.yaml pins its commit), per
# Ribo-seq library: P-site offsets from the reads on the human transcripts,
# then reads per human CDS, and per spike-in (yeast) CDS with the human
# offsets. Inputs: align.smk's split human / spike-in BAMs and the
# experiment's RNA-seq-filtered GTF, and the genome fastas and spike-in GTF of
# config references:. Only reads of config read_lengths count, the same window
# in every library. In results/ribokit/<exp>/, per library and species (human,
# yeast):
#   <sample>.<species>.quant.tsv     reads and ritpm per CDS (Name Length
#                                    EffectiveLength ritpm NumReads)
#   <sample>.<species>.offsets.tsv   P-site offsets per read length and phase
#                                    (yeast: the human ones it used)
#   <sample>.<species>.ties.tsv      transcripts with identical CDSs, grouped
#   <sample>.<species>.stats.tsv     annotation, read and fit counts
# Later steps read the quants (RIBOKIT_QUANT). ribokit fails by itself if its
# fit does not converge.

import re

_rl = config["read_lengths"]
if not (isinstance(_rl, list) and len(_rl) == 2 and all(isinstance(x, int) for x in _rl) and _rl[0] <= _rl[1]):
    raise ValueError(f"config read_lengths: needs [lo, hi]; got {_rl}")

RIBOKIT_DIR = f"{RESULTS_DIR}/ribokit"
RIBOKIT_OUT = f"{RIBOKIT_DIR}/{{exp}}/{{sample}}"
RIBOKIT_QUANT = f"{RIBOKIT_OUT}.{{species}}.quant.tsv"
RIBOKIT_CONSTRAINTS = dict(
    exp="|".join(map(re.escape, EXPERIMENT_NAMES)),
    sample="|".join(map(re.escape, RIBO_SAMPLES)),
)
RIBOKIT_LOG = f"{LOG_DIR}/ribokit/{{exp}}/{{sample}}"
RIBOKIT_CMD = f"ribokit quant --read-lengths {_rl[0]}-{_rl[1]} "


rule ribokit_human_morfs:
    input:
        bam=f"{ALIGN_DIR}/split_bam/transcriptome/human/{{sample}}.bam",
        fasta=config["references"]["human_genome_fa"],
        gtf=f"{ALIGN_DIR}/reference/human_transcriptome.rnaseq_filtered.gtf",
    output:
        quant=f"{RIBOKIT_OUT}.human.quant.tsv",
        offsets=f"{RIBOKIT_OUT}.human.offsets.tsv",
        ties=f"{RIBOKIT_OUT}.human.ties.tsv",
        stats=f"{RIBOKIT_OUT}.human.stats.tsv",
    wildcard_constraints:
        **RIBOKIT_CONSTRAINTS,
    params:
        prefix=f"{RIBOKIT_OUT}.human",
    log:
        f"{RIBOKIT_LOG}.human_morfs.log",
    conda:
        "../envs/ribokit.yaml"
    shell:
        RIBOKIT_CMD + "--bam {input.bam} --gtf {input.gtf} --fasta {input.fasta} "
        "--out-prefix {params.prefix} > {log} 2>&1"


rule ribokit_yeast_morfs:
    input:
        bam=f"{ALIGN_DIR}/split_bam/transcriptome/spike_in/{{sample}}.bam",
        fasta=config["references"]["yeast_genome_fa"],
        gtf=config["references"]["yeast_gtf"],
        offsets=rules.ribokit_human_morfs.output.offsets,
    output:
        quant=f"{RIBOKIT_OUT}.yeast.quant.tsv",
        offsets=f"{RIBOKIT_OUT}.yeast.offsets.tsv",
        ties=f"{RIBOKIT_OUT}.yeast.ties.tsv",
        stats=f"{RIBOKIT_OUT}.yeast.stats.tsv",
    wildcard_constraints:
        **RIBOKIT_CONSTRAINTS,
    params:
        prefix=f"{RIBOKIT_OUT}.yeast",
    log:
        f"{RIBOKIT_LOG}.yeast_morfs.log",
    conda:
        "../envs/ribokit.yaml"
    shell:
        RIBOKIT_CMD + "--bam {input.bam} --gtf {input.gtf} --fasta {input.fasta} --offsets {input.offsets} "
        "--out-prefix {params.prefix} > {log} 2>&1"
