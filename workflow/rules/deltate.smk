# The deltaTE model (Chotani et al. 2019), per experiment: DESeq2 on its
# Ribo-seq and RNA-seq libraries together, the change in translational
# efficiency being the interaction of condition (config conditions: reference
# -> treatment) and assay (RNA-seq -> Ribo-seq). Inputs: ribokit.smk's quants
# (RIBOKIT_QUANT: human for the model, yeast for the spike-in size factors) and
# rna_correction.smk's (RNA_CORRECTION_QUANT). Each library's assay and
# condition come from samples.csv. The scripts read and write
# results/deltaTE/<exp>/:
#   sizeFactors_yeast_ribo.rds       calc_sizefactors.R: Ribo-seq size factors
#                                    from the spike-in
#   sampleTable.*, txi*              import_txAbundance.R: the libraries in one
#                                    tximport object
#   deseq_res_deltaTE_<norm>.*       run_deltaTE.R: the change in TE
#   deseq_res_diff{ribo,total}_<norm>.*   run_deseq.R: each assay alone
#   plots/MAplot_diff_*.pdf
# <norm>: yeastnorm, the Ribo-seq size factors from the spike-in; autonorm,
# estimated within the Ribo-seq libraries (RNA-seq: autonorm only).

import re

# [reference, treatment]
DELTATE_CONDITIONS = [config["conditions"]["reference"], config["conditions"]["treatment"]]


def _check_samples(samples):
    """Every experiment run has exactly the two conditions in both assays. The
    scripts take each library's assay and condition from samples.csv, never
    from its name (wf-riboseq-align checks read_type and sample_id)."""
    errors = []
    for exp, sub in samples[samples["experiment"].isin(EXPERIMENT_NAMES)].groupby("experiment"):
        for read_type in ("ribo", "rna"):
            conds = set(sub.loc[sub["read_type"] == read_type, "condition"])
            if conds != set(DELTATE_CONDITIONS):
                errors.append(f"{exp}: {read_type} has conditions {sorted(conds)}, needs {DELTATE_CONDITIONS}")
    if errors:
        raise ValueError(f"{config['samples']}:\n  " + "\n  ".join(errors))


_check_samples(SAMPLES)


DELTATE_DIR = f"{RESULTS_DIR}/deltaTE"
DELTATE_OUT = f"{DELTATE_DIR}/{{exp}}"
DELTATE_CONSTRAINTS = dict(exp="|".join(map(re.escape, EXPERIMENT_NAMES)))
DELTATE_LOG = f"{LOG_DIR}/deltaTE/{{exp}}"


def deltate_quants(path, read_type, **wildcards):
    """An experiment's quants of one read type, in samples.csv's order:
    calc_sizefactors.R scales every Ribo-seq library to the first."""
    return lambda wc: expand(path, exp=wc.exp, sample=samples_of(wc.exp, read_type), **wildcards)


def deltate_results(kind, norms):
    return [f"{DELTATE_OUT}/deseq_res_{kind}_{norm}.{ext}" for norm in norms for ext in ("tsv", "rds")]


def deltate_plots(kind, norms):
    return [f"{DELTATE_OUT}/plots/MAplot_diff_{kind}_{norm}.pdf" for norm in norms]


rule deltate_calc_sizefactors:
    input:
        yeast=deltate_quants(RIBOKIT_QUANT, "ribo", species="yeast"),
        script=workflow.source_path("../scripts/calc_sizefactors.R"),
    output:
        f"{DELTATE_OUT}/sizeFactors_yeast_ribo.rds",
    wildcard_constraints:
        **DELTATE_CONSTRAINTS,
    params:
        dir=DELTATE_OUT,
        samples=lambda wc: " ".join([str(len(samples_of(wc.exp, "ribo")))] + samples_of(wc.exp, "ribo")),
    log:
        f"{DELTATE_LOG}.calc_sizefactors.log",
    conda:
        "../envs/deseq2.yaml"
    shell:
        "Rscript {input.script} {params.dir} {params.samples} {input.yeast} > {log} 2>&1"


rule deltate_import_txabundance:
    input:
        ribo=deltate_quants(RIBOKIT_QUANT, "ribo", species="human"),
        rna=deltate_quants(RNA_CORRECTION_QUANT, "rna"),
        script=workflow.source_path("../scripts/import_txAbundance.R"),
    output:
        sample_table=f"{DELTATE_OUT}/sampleTable.rds",
        sample_table_tsv=f"{DELTATE_OUT}/sampleTable.tsv",
        txi=f"{DELTATE_OUT}/txi.rds",
        txi_tsv=[f"{DELTATE_OUT}/txi_{m}.tsv" for m in ("counts", "abundance", "length")],
    wildcard_constraints:
        **DELTATE_CONSTRAINTS,
    params:
        dir=DELTATE_OUT,
        conditions=" ".join(DELTATE_CONDITIONS),
        samples=lambda wc: " ".join([str(len(samples_of(wc.exp, "ribo"))), str(len(samples_of(wc.exp, "rna")))]
                                    + samples_of(wc.exp, "ribo") + samples_of(wc.exp, "rna")),
        sample_conditions=lambda wc: " ".join(SAMPLES.loc[samples_of(wc.exp, "ribo") + samples_of(wc.exp, "rna"),
                                                          "condition"]),
    log:
        f"{DELTATE_LOG}.import_txabundance.log",
    conda:
        "../envs/deseq2.yaml"
    shell:
        "Rscript {input.script} {params.dir} {params.conditions} {params.samples} {params.sample_conditions} "
        "{input.ribo} {input.rna} > {log} 2>&1"


rule deltate_run_deltaTE:
    input:
        sample_table=rules.deltate_import_txabundance.output.sample_table,
        txi=rules.deltate_import_txabundance.output.txi,
        size_factors=rules.deltate_calc_sizefactors.output[0],
        script=workflow.source_path("../scripts/run_deltaTE.R"),
    output:
        results=deltate_results("deltaTE", ("yeastnorm", "autonorm")),
        size_factors=[f"{DELTATE_OUT}/sizeFactors_human_{x}.rds" for x in ("total", "ribo_auto")],
        plots=deltate_plots("deltaTE", ("yeastnorm", "autonorm")),
    wildcard_constraints:
        **DELTATE_CONSTRAINTS,
    params:
        dir=DELTATE_OUT,
    log:
        f"{DELTATE_LOG}.run_deltaTE.log",
    conda:
        "../envs/deseq2.yaml"
    shell:
        "Rscript {input.script} {params.dir} > {log} 2>&1"


rule deltate_run_deseq_total:
    input:
        sample_table=rules.deltate_import_txabundance.output.sample_table,
        txi=rules.deltate_import_txabundance.output.txi,
        script=workflow.source_path("../scripts/run_deseq.R"),
    output:
        results=deltate_results("difftotal", ("autonorm",)),
        plots=deltate_plots("total", ("autonorm",)),
    wildcard_constraints:
        **DELTATE_CONSTRAINTS,
    params:
        dir=DELTATE_OUT,
    log:
        f"{DELTATE_LOG}.run_deseq_total.log",
    conda:
        "../envs/deseq2.yaml"
    shell:
        "Rscript {input.script} {params.dir} total > {log} 2>&1"


rule deltate_run_deseq_ribo:
    input:
        sample_table=rules.deltate_import_txabundance.output.sample_table,
        txi=rules.deltate_import_txabundance.output.txi,
        size_factors=rules.deltate_calc_sizefactors.output[0],
        script=workflow.source_path("../scripts/run_deseq.R"),
    output:
        results=deltate_results("diffribo", ("autonorm", "yeastnorm")),
        plots=deltate_plots("ribo", ("autonorm", "yeastnorm")),
    wildcard_constraints:
        **DELTATE_CONSTRAINTS,
    params:
        dir=DELTATE_OUT,
    log:
        f"{DELTATE_LOG}.run_deseq_ribo.log",
    conda:
        "../envs/deseq2.yaml"
    shell:
        "Rscript {input.script} {params.dir} ribo > {log} 2>&1"
