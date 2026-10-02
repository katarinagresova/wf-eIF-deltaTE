# The RNA-seq library correction (config rna_correction:), per experiment: a
# library whose counts lean with transcript length and GC against its
# replicates is corrected (scripts/rna_correction.R says how). Input: align.smk's
# salmon/<sample>/quant.rnaseq_filtered.sf (salmon's quant.sf without the
# transcripts the RNA-seq filter removed). In results/rna_correction/<exp>/:
#   <sample>_quant.sf   the corrected libraries rewritten, every other file copied
#   scan.tsv            every fit, round by round; factors.tsv
# Later steps read the quants (RNA_CORRECTION_QUANT). Only the RNA-seq is
# touched: nothing computed from the Ribo-seq alone depends on it.

RNA_CORRECTION = config["rna_correction"]
RNA_CORRECTION_DIR = f"{RESULTS_DIR}/rna_correction"
RNA_CORRECTION_FEATURES = f"{RNA_CORRECTION_DIR}/transcript_features.tsv"
RNA_CORRECTION_QUANT = f"{RNA_CORRECTION_DIR}/{{exp}}/{{sample}}_quant.sf"


# length and GC of every transcript of the reference salmon quantifies against
rule rna_correction_features:
    input:
        fa=config["align"]["human_transcriptome_fa"],
        script=workflow.source_path("../scripts/transcript_features.py"),
    output:
        RNA_CORRECTION_FEATURES,
    log:
        f"{LOG_DIR}/rna_correction/transcript_features.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.fa} {output} 2> {log}"


# One rule per experiment (`name:` in a loop): its outputs are a list of that
# experiment's libraries.
for _exp in EXPERIMENT_NAMES:
    _dir = f"{RNA_CORRECTION_DIR}/{_exp}"
    _rna = samples_of(_exp, "rna")

    rule:
        name:
            f"rna_correction_{_exp}"
        input:
            quant=[f"{ALIGN_DIR.format(exp=_exp)}/salmon/{s}/quant.rnaseq_filtered.sf" for s in _rna],
            features=RNA_CORRECTION_FEATURES,
            script=workflow.source_path("../scripts/rna_correction.R"),
        output:
            quant=[f"{_dir}/{s}_quant.sf" for s in _rna],
            scan=f"{_dir}/scan.tsv",
            factors=f"{_dir}/factors.tsv",
        params:
            settings=f"{RNA_CORRECTION['max_sd']} {RNA_CORRECTION['min_mean_count']} {RNA_CORRECTION['spline_df']}",
            n=len(_rna),
            conditions=" ".join(SAMPLES.loc[_rna, "condition"]),
        log:
            f"{LOG_DIR}/rna_correction/{_exp}.log",
        conda:
            "../envs/rna_correction.yaml"
        shell:
            "Rscript {input.script} {input.features} {params.settings} {output.scan} {output.factors} {params.n} "
            "{input.quant} {output.quant} {params.conditions} > {log} 2>&1"
