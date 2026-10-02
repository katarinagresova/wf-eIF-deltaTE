# The eIF master table (workflow/scripts/master_table.py), the table the
# paper reads through wf-transmod-data (its HCT116_ISOFORMS_TE: logTE). Same
# columns and row order as Frederick's results/20260324-eif-master.csv, built
# from this workflow's deltaTE results and quants instead of his folders.
#
# Validated on his own 10 folders: same 158,284 rows, columns and row order as
# 20260324-eif-master.csv, every column equal except log2FC_TE / padj_TE, by at
# most 4e-5 / 3e-4 with identical hits in all 10 experiments. That is his own
# 2026-07-16 rerun of the deltaTE tables (his folders now hold those), not the port.

# His notebook's experiment order (factors, then 4h before 8h); any other target after them
MASTER_TARGET_ORDER = ["eIF4G1", "eIF4G2", "eIF4G3", "eIF4E", "eIF3d"]


def master_experiments():
    def key(exp):
        target, timepoint = exp.rsplit("_", 1)
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
        gtf=config["align"]["human_transcriptome_gtf"],
        script=workflow.source_path("../scripts/master_table.py"),
    output:
        f"{RESULTS_DIR}/master/eif-master.csv",
    params:
        args=" ".join(" ".join(["--experiment", e, *master_tables(e)]
                               + [x for lib in master_libraries(e) for x in ("--library", e, *lib)])
                      for e in master_experiments()),
    log:
        f"{LOG_DIR}/master.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.gtf} {output} {params.args} > {log} 2>&1"
