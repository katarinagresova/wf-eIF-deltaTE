# The matched total RNA-seq (paired-end), per experiment: adapter trimming
# (config rnaseq_trim:), salmon on the whole
# (collapsed) human transcriptome, the genome as decoys, then the expression
# filter (config autofilter:), whose blacklist is applied to what is counted:
# the GTF ribokit reads and the quants the RNA-seq library correction reads.
# The Ribo-seq alignment does not use the blacklist (align.smk aligns to the
# whole transcriptome); filter_ribo_bam writes a copy of the human BAMs
# without the blacklisted transcripts, for tools other than this workflow's.
# The salmon index does not depend on the experiment:
# results/rnaseq/salmon_index/. In results/rnaseq/<exp>/:
#   salmon/<sample>/                         salmon's output (quant.sf, its logs)
#   salmon/<sample>/quant.rnaseq_filtered.sf quant.sf minus the blacklist (RNASEQ_QUANT)
#   rnaseq_filter_blacklist_txid.txt         the transcripts below the floor
#   human_transcriptome.rnaseq_filtered.gtf  the collapsed GTF minus them (RNASEQ_GTF)
#   ribo_bam/<ribo sample>.human.rnaseq_filtered.bam(.bai)
#                                            align.smk's human Ribo-seq BAM minus them (RNASEQ_RIBO_BAM)
#   trimmed/<sample>_R{1,2}.fastq.gz         the adapter-trimmed fastqs (temp)
#   fastqc/<sample>_R{1,2}_fastqc.{html,zip} FastQC of the raw fastqs
#   fastqc_trimmed/<sample>_R{1,2}_fastqc.{html,zip}
#                                            FastQC of the trimmed fastqs
#   multiqc/multiqc_report.html              MultiQC: cutadapt, salmon, both FastQCs

import re

RNASEQ_DIR = f"{RESULTS_DIR}/rnaseq"
RNASEQ_LOG = f"{LOG_DIR}/rnaseq"
SALMON_INDEX_DIR = f"{RNASEQ_DIR}/salmon_index"
RNASEQ_BLACKLIST = f"{RNASEQ_DIR}/{{exp}}/rnaseq_filter_blacklist_txid.txt"
RNASEQ_GTF = f"{RNASEQ_DIR}/{{exp}}/human_transcriptome.rnaseq_filtered.gtf"
RNASEQ_QUANT = f"{RNASEQ_DIR}/{{exp}}/salmon/{{sample}}/quant.rnaseq_filtered.sf"
RNASEQ_RIBO_BAM = f"{RNASEQ_DIR}/{{exp}}/ribo_bam/{{sample}}.human.rnaseq_filtered.bam"
RNASEQ_CONSTRAINTS = dict(
    exp=EXPERIMENT_RE,
    sample="|".join(map(re.escape, RNA_SAMPLES)),
)


def _rnaseq_files(pattern):
    """`pattern` ({exp}, {sample}) for every RNA-seq library of the job's experiment."""
    return lambda wc: [pattern.format(exp=wc.exp, sample=s) for s in samples_of(wc.exp, "rna")]


# salmon's decoy-aware reference: the transcripts, then every genome sequence
# (decoys must come last), and the decoys' names.
rule salmon_gentrome:
    input:
        transcriptome=COLLAPSED_FA,
        genome=config["references"]["human_genome_fa"],
    output:
        gentrome=temp(f"{SALMON_INDEX_DIR}/gentrome.fa"),
        decoys=f"{SALMON_INDEX_DIR}/decoys.txt",
    log:
        f"{RNASEQ_LOG}/salmon_index/gentrome.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        r"""
        exec 2> {log}
        grep '^>' {input.genome} | cut -c2- | cut -d' ' -f1 > {output.decoys}
        cat {input.transcriptome} {input.genome} > {output.gentrome}
        """


# A fragment that aligns better to the genome than to any transcript (intronic,
# intergenic, an unannotated copy) is counted for no transcript. By default
# salmon drops a transcript whose sequence duplicates another's from the index
# and quant.sf, and filter_rnaseq would blacklist it as TPM 0:
# --keepDuplicates.
rule salmon_index:
    input:
        gentrome=f"{SALMON_INDEX_DIR}/gentrome.fa",
        decoys=f"{SALMON_INDEX_DIR}/decoys.txt",
    output:
        directory(f"{SALMON_INDEX_DIR}/transcriptome_genome_decoys"),
    log:
        f"{RNASEQ_LOG}/salmon_index/index.log",
    conda:
        "../envs/salmon.yaml"
    threads: 16
    resources:
        mem_mb=24000,
    shell:
        "salmon index -p {threads} -t {input.gentrome} -d {input.decoys} -k 31 --keepDuplicates "
        "-i {output} > {log} 2>&1"


# eIF3d's fastqs come untrimmed (1-5% of reads reach the adapter), the other
# targets' trimmed. Every library is trimmed the same way here.
# envs/cutadapt.yaml is wf-riboseq-align's at ALIGN_COMMIT byte for byte, so
# the align module's built env is reused.
rule cutadapt_rna:
    input:
        fastq1=lambda wc: SAMPLES.loc[wc.sample, "fastq_1"],
        fastq2=lambda wc: SAMPLES.loc[wc.sample, "fastq_2"],
    output:
        fastq1=temp(f"{RNASEQ_DIR}/{{exp}}/trimmed/{{sample}}_R1.fastq.gz"),
        fastq2=temp(f"{RNASEQ_DIR}/{{exp}}/trimmed/{{sample}}_R2.fastq.gz"),
    wildcard_constraints:
        **RNASEQ_CONSTRAINTS,
    params:
        adapter=config["rnaseq_trim"]["adapter"],
        minlength=config["rnaseq_trim"]["minlength"],
        qualcutoff=config["rnaseq_trim"]["qualcutoff"],
    log:
        f"{RNASEQ_LOG}/{{exp}}/cutadapt/{{sample}}.log",
    conda:
        "../envs/cutadapt.yaml"
    threads: 4
    shell:
        # -O 6 (as wf-riboseq-align's cutadapt_reads): default -O 3 lets a 3 nt
        # chance match to the adapter trim real insert bases off the 3' end;
        # -O 6 is long enough to be a real adapter. -O and the quality cutoff
        # apply to both mates.
        "cutadapt --cores {threads} -a {params.adapter} -A {params.adapter} -O 6 "
        "--minimum-length {params.minlength} --quality-cutoff {params.qualcutoff} "
        "-o {output.fastq1} -p {output.fastq2} "
        "{input.fastq1} {input.fastq2} > {log} 2>&1"


rule salmon_quant:
    input:
        fastq1=f"{RNASEQ_DIR}/{{exp}}/trimmed/{{sample}}_R1.fastq.gz",
        fastq2=f"{RNASEQ_DIR}/{{exp}}/trimmed/{{sample}}_R2.fastq.gz",
        index=f"{SALMON_INDEX_DIR}/transcriptome_genome_decoys",
    output:
        f"{RNASEQ_DIR}/{{exp}}/salmon/{{sample}}/quant.sf",
    wildcard_constraints:
        **RNASEQ_CONSTRAINTS,
    params:
        out_dir=f"{RNASEQ_DIR}/{{exp}}/salmon/{{sample}}",
    log:
        f"{RNASEQ_LOG}/{{exp}}/salmon/{{sample}}.log",
    conda:
        "../envs/salmon.yaml"
    threads: 8
    resources:
        mem_mb=24000,
    shell:
        "salmon quant -p {threads} -i {input.index} -l A -1 {input.fastq1} -2 {input.fastq2} "
        "--seqBias --gcBias --output {params.out_dir} > {log} 2>&1"


rule filter_rnaseq:
    input:
        fa=COLLAPSED_FA,
        quant=_rnaseq_files(f"{RNASEQ_DIR}/{{exp}}/salmon/{{sample}}/quant.sf"),
        script=workflow.source_path("../scripts/autofilter_rnaseq.py"),
    output:
        RNASEQ_BLACKLIST,
    wildcard_constraints:
        exp=RNASEQ_CONSTRAINTS["exp"],
    params:
        min_tpm=config["autofilter"]["min_tpm"],
        min_samples=config["autofilter"]["min_samples"],
    log:
        f"{RNASEQ_LOG}/{{exp}}/filter_rnaseq.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.fa} {output} "
        "{params.min_tpm} {params.min_samples} {input.quant} 2> {log}"


rule filter_gtf:
    input:
        gtf=COLLAPSED_GTF,
        blacklist=RNASEQ_BLACKLIST,
        script=workflow.source_path("../scripts/filter_gtf.py"),
    output:
        RNASEQ_GTF,
    wildcard_constraints:
        exp=RNASEQ_CONSTRAINTS["exp"],
    log:
        f"{RNASEQ_LOG}/{{exp}}/filter_gtf.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.gtf} {input.blacklist} {output} 2> {log}"


# The RNA-seq quantification on the same transcripts: salmon's quant.sf minus
# the blacklisted rows (filter_quant.py says what is and is not changed).
rule filter_quant:
    input:
        quant=f"{RNASEQ_DIR}/{{exp}}/salmon/{{sample}}/quant.sf",
        blacklist=RNASEQ_BLACKLIST,
        script=workflow.source_path("../scripts/filter_quant.py"),
    output:
        RNASEQ_QUANT,
    wildcard_constraints:
        **RNASEQ_CONSTRAINTS,
    log:
        f"{RNASEQ_LOG}/{{exp}}/filter_quant/{{sample}}.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.quant} {input.blacklist} {output} 2> {log}"


# The Ribo-seq alignments on the same transcripts, for tools other than this
# workflow's (ribokit quant already drops what is not in RNASEQ_GTF): the human
# split BAM minus the blacklisted transcripts (filter_bam.py says what is and
# is not changed).
rule filter_ribo_bam:
    input:
        bam=f"{ALIGN_DIR}/split_bam/transcriptome/human/{{sample}}.bam",
        blacklist=RNASEQ_BLACKLIST,
        script=workflow.source_path("../scripts/filter_bam.py"),
    output:
        bam=RNASEQ_RIBO_BAM,
        bai=f"{RNASEQ_RIBO_BAM}.bai",
    wildcard_constraints:
        exp=RNASEQ_CONSTRAINTS["exp"],
        sample="|".join(map(re.escape, RIBO_SAMPLES)),
    log:
        f"{RNASEQ_LOG}/{{exp}}/filter_ribo_bam/{{sample}}.log",
    conda:
        "../envs/pysam.yaml"
    shell:
        "python {input.script} {input.bam} {input.blacklist} {output.bam} 2> {log}"


# FastQC of the raw RNA-seq, R1 and R2. Reports only; no rule but MultiQC reads
# them. FastQC names its report after the input file, so each job links its
# fastq into its own tmp dir as <name>.fastq.gz and moves the report out. The
# tmp dir also takes FastQC's temporary files (--dir): compute nodes' /tmp can
# be full. (wf-riboseq-align's FASTQC_SHELL.)
rule fastqc_rna_raw:
    input:
        fastq=lambda wc: SAMPLES.loc[wc.sample, f"fastq_{wc.mate}"],
    output:
        html=f"{RNASEQ_DIR}/{{exp}}/fastqc/{{sample}}_R{{mate}}_fastqc.html",
        zip=f"{RNASEQ_DIR}/{{exp}}/fastqc/{{sample}}_R{{mate}}_fastqc.zip",
    wildcard_constraints:
        mate="[12]",
        **RNASEQ_CONSTRAINTS,
    params:
        name="{sample}_R{mate}",
        tmp=f"{RNASEQ_DIR}/{{exp}}/fastqc/{{sample}}_R{{mate}}_tmp",
    log:
        f"{RNASEQ_LOG}/{{exp}}/fastqc/{{sample}}_R{{mate}}.log",
    conda:
        "../envs/fastqc.yaml"
    shell:
        "rm -rf {params.tmp} && mkdir -p {params.tmp} && "
        "ln -s $(readlink -f {input.fastq}) {params.tmp}/{params.name}.fastq.gz && "
        "fastqc --quiet --threads 1 --dir {params.tmp} --outdir {params.tmp} "
        "{params.tmp}/{params.name}.fastq.gz > {log} 2>&1 && "
        "mv {params.tmp}/{params.name}_fastqc.html {output.html} && "
        "mv {params.tmp}/{params.name}_fastqc.zip {output.zip} && "
        "rm -rf {params.tmp}"


# FastQC of the trimmed RNA-seq, what salmon reads. The trimmed fastqs are
# already named <sample>_R<mate>.fastq.gz, so no link; the tmp dir is
# FastQC's --dir only.
rule fastqc_rna_trimmed:
    input:
        fastq=f"{RNASEQ_DIR}/{{exp}}/trimmed/{{sample}}_R{{mate}}.fastq.gz",
    output:
        html=f"{RNASEQ_DIR}/{{exp}}/fastqc_trimmed/{{sample}}_R{{mate}}_fastqc.html",
        zip=f"{RNASEQ_DIR}/{{exp}}/fastqc_trimmed/{{sample}}_R{{mate}}_fastqc.zip",
    wildcard_constraints:
        mate="[12]",
        **RNASEQ_CONSTRAINTS,
    params:
        tmp=f"{RNASEQ_DIR}/{{exp}}/fastqc_trimmed/{{sample}}_R{{mate}}_tmp",
    log:
        f"{RNASEQ_LOG}/{{exp}}/fastqc_trimmed/{{sample}}_R{{mate}}.log",
    conda:
        "../envs/fastqc.yaml"
    shell:
        "rm -rf {params.tmp} && mkdir -p {params.tmp} && "
        "fastqc --quiet --threads 1 --dir {params.tmp} --outdir $(dirname {output.html}) "
        "{input.fastq} > {log} 2>&1 && "
        "rm -rf {params.tmp}"


# One MultiQC report per experiment: cutadapt, salmon and the RNA-seq FastQC,
# raw and trimmed (sections and sample names: workflow/multiqc_rnaseq.yaml).
# The Ribo-seq's is align.smk's (results/align/<exp>/multiqc/).
rule multiqc_rnaseq:
    input:
        fastqc=lambda wc: [f"{RNASEQ_DIR}/{wc.exp}/{d}/{s}_R{mate}_fastqc.zip"
                           for d in ("fastqc", "fastqc_trimmed")
                           for s in samples_of(wc.exp, "rna") for mate in (1, 2)],
        cutadapt=_rnaseq_files(f"{RNASEQ_LOG}/{{exp}}/cutadapt/{{sample}}.log"),
        salmon=_rnaseq_files(f"{RNASEQ_DIR}/{{exp}}/salmon/{{sample}}/quant.sf"),
        config=workflow.source_path("../multiqc_rnaseq.yaml"),
    output:
        html=f"{RNASEQ_DIR}/{{exp}}/multiqc/multiqc_report.html",
        data=directory(f"{RNASEQ_DIR}/{{exp}}/multiqc/multiqc_data"),
    wildcard_constraints:
        exp=RNASEQ_CONSTRAINTS["exp"],
    # salmon's report is its output directory, not quant.sf
    params:
        salmon_dirs=_rnaseq_files(f"{RNASEQ_DIR}/{{exp}}/salmon/{{sample}}"),
    log:
        f"{RNASEQ_LOG}/{{exp}}/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc --force --config {input.config} --outdir $(dirname {output.html}) "
        "{input.fastqc} {input.cutadapt} {params.salmon_dirs} > {log} 2>&1"
