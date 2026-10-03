# The deltaTE model (Chothani et al. 2019), per experiment: DESeq2 on its
# Ribo-seq and RNA-seq libraries together, the change in translational
# efficiency being the interaction of condition (config conditions: reference
# -> treatment) and assay (RNA-seq -> Ribo-seq); each assay is also fitted
# alone. Inputs: ribokit.smk's quants (RIBOKIT_QUANT: human for the model,
# yeast for the spike-in size factors) and rna_correction.smk's
# (RNA_CORRECTION_QUANT). Each library's assay and condition come from
# samples.csv. In results/deltaTE/<exp>/:
#   sampleTable.tsv, counts.tsv, length.tsv, spike_in_reads.tsv
#                                deltate_inputs.py: the libraries as DESeq2's
#                                input, the spike-in reads per Ribo-seq library
#   deseq_res_<model>.tsv        deseq2.R: deltaTE_<norm> (the change in TE),
#   plots/MAplot_diff_*.pdf      diffribo_<norm>, difftotal_autonorm
# <norm>: yeastnorm, the Ribo-seq size factors from the spike-in; autonorm,
# estimated within the Ribo-seq libraries (RNA-seq: autonorm only).

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
DELTATE_LOG = f"{LOG_DIR}/deltaTE/{{exp}}"
# deseq2.R's models: deseq_res_<model>.tsv, plots/MAplot_diff_<plot>.pdf
DELTATE_MODELS = {"deltaTE_yeastnorm": "deltaTE_yeastnorm", "deltaTE_autonorm": "deltaTE_autonorm",
                  "diffribo_yeastnorm": "ribo_yeastnorm", "diffribo_autonorm": "ribo_autonorm",
                  "difftotal_autonorm": "total_autonorm"}


def deltate_quants(path, read_type, **wildcards):
    """An experiment's quants of one read type, in samples.csv's order:
    deseq2.R scales every Ribo-seq library's spike-in to the first's."""
    return lambda wc: expand(path, exp=wc.exp, sample=samples_of(wc.exp, read_type), **wildcards)


def deltate_libraries(wc, input):
    """deltate_inputs.py's --ribo / --rna arguments: per library its sample id,
    condition and quant(s)."""
    args = [x for s, human, yeast in zip(samples_of(wc.exp, "ribo"), input.human, input.yeast)
            for x in ("--ribo", s, SAMPLES.at[s, "condition"], human, yeast)]
    args += [x for s, quant in zip(samples_of(wc.exp, "rna"), input.rna)
             for x in ("--rna", s, SAMPLES.at[s, "condition"], quant)]
    return " ".join(args)


rule deltate_inputs:
    input:
        human=deltate_quants(RIBOKIT_QUANT, "ribo", species="human"),
        yeast=deltate_quants(RIBOKIT_QUANT, "ribo", species="yeast"),
        rna=deltate_quants(RNA_CORRECTION_QUANT, "rna"),
        script=workflow.source_path("../scripts/deltate_inputs.py"),
    output:
        [f"{DELTATE_OUT}/{x}.tsv" for x in ("sampleTable", "counts", "length", "spike_in_reads")],
    wildcard_constraints:
        exp=EXPERIMENT_RE,
    params:
        dir=DELTATE_OUT,
        libraries=deltate_libraries,
    log:
        f"{DELTATE_LOG}.inputs.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {params.dir} {params.libraries} > {log} 2>&1"


rule deltate_deseq2:
    input:
        tables=rules.deltate_inputs.output,
        script=workflow.source_path("../scripts/deseq2.R"),
    output:
        results=[f"{DELTATE_OUT}/deseq_res_{model}.tsv" for model in DELTATE_MODELS],
        plots=[f"{DELTATE_OUT}/plots/MAplot_diff_{plot}.pdf" for plot in DELTATE_MODELS.values()],
    wildcard_constraints:
        exp=EXPERIMENT_RE,
    params:
        dir=DELTATE_OUT,
        conditions=" ".join(DELTATE_CONDITIONS),
    log:
        f"{DELTATE_LOG}.deseq2.log",
    conda:
        "../envs/deseq2.yaml"
    shell:
        "Rscript {input.script} {params.dir} {params.conditions} > {log} 2>&1"
