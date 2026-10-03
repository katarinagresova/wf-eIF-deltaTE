# The eIF master table (workflow/scripts/master_table.py), the table the
# paper reads through wf-transmod-data (its HCT116_ISOFORMS_TE: logTE), from
# the deltaTE results and quants of every experiment run.

# The experiments' order: these targets first, in this order, any other target
# after them by name; within a target by timepoint
MASTER_TARGET_ORDER = ["eIF4G1", "eIF4G2", "eIF4G3", "eIF4E", "eIF3d"]


def master_experiments():
    def key(exp):
        target, timepoint = SAMPLES.loc[SAMPLES["experiment"] == exp, ["target", "timepoint"]].iloc[0]
        rank = MASTER_TARGET_ORDER.index(target) if target in MASTER_TARGET_ORDER else len(MASTER_TARGET_ORDER)
        return rank, target, timepoint
    return sorted(EXPERIMENT_NAMES, key=key)


def master_tables(exp):
    """TE, RPF and RNA: the experiment's DESeq2 tables the master table takes."""
    return [f"{DELTATE_DIR}/{exp}/deseq_res_{k}.tsv"
            for k in ("deltaTE_yeastnorm", "diffribo_yeastnorm", "difftotal_autonorm")]


def master_libraries(exp):
    """Per library of the experiment: read_type, condition, rep (samples.csv) and its quant."""
    quants = [(s, RIBOKIT_QUANT.format(exp=exp, sample=s, species="human")) for s in samples_of(exp, "ribo")]
    quants += [(s, RNA_CORRECTION_QUANT.format(exp=exp, sample=s)) for s in samples_of(exp, "rna")]
    return [(SAMPLES.at[s, "read_type"], SAMPLES.at[s, "condition"], SAMPLES.at[s, "rep"], q) for s, q in quants]


rule master_table:
    input:
        tables=[t for e in master_experiments() for t in master_tables(e)],
        quants=[lib[3] for e in master_experiments() for lib in master_libraries(e)],
        gtf=COLLAPSED_GTF,
        script=workflow.source_path("../scripts/master_table.py"),
    output:
        f"{RESULTS_DIR}/master/eif-master.csv",
    params:
        args=" ".join(["--conditions", *DELTATE_CONDITIONS]
                      + [x for e in master_experiments() for x in ["--experiment", e, *master_tables(e)]]
                      + [x for e in master_experiments() for lib in master_libraries(e)
                         for x in ("--library", e, *lib)]),
    log:
        f"{LOG_DIR}/master.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.gtf} {output} {params.args} > {log} 2>&1"


# The master table restricted to the transcripts expressed in the reference
# condition (config master_filter, workflow/scripts/filter_master.py): per
# experiment, Ribo-seq ritpm and RNA-seq TPM above their floors in every
# reference library. After DESeq2, so padj is the full table's.
rule master_table_filtered:
    input:
        master=rules.master_table.output[0],
        script=workflow.source_path("../scripts/filter_master.py"),
    output:
        f"{RESULTS_DIR}/master/eif-master.filtered.csv",
    params:
        args=" ".join(map(str, [DELTATE_CONDITIONS[0], config["master_filter"]["min_ribo_tpm"],
                                config["master_filter"]["min_rna_tpm"]])),
    log:
        f"{LOG_DIR}/master_filtered.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.master} {output} {params.args} > {log} 2>&1"
