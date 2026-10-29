# The Ribo-seq QC panel, per experiment: from ribokit.smk's <sample>.human.qc.tsv
# (reads per read length, P-site region and 5' end phase), the footprint length
# distribution, the P-site region and the 5' end phase in the CDS, per read
# length, in the config read_lengths window. In results/rpf_qc/<exp>/:
#   <exp>_rpf_qc.pdf     the three plots in a row, a line or dot per library

RPF_QC_DIR = f"{RESULTS_DIR}/rpf_qc"
RPF_QC_OUT = f"{RPF_QC_DIR}/{{exp}}"
RPF_QC_LOG = f"{LOG_DIR}/rpf_qc/{{exp}}"


def rpf_qc_tables(wc):
    """An experiment's ribokit qc tables, in samples.csv's order."""
    return expand(f"{RIBOKIT_OUT}.human.qc.tsv", exp=wc.exp, sample=samples_of(wc.exp, "ribo"))


def rpf_qc_libraries(wc, input):
    """plot_rpf_qc.py's --lib arguments: per library its sample id, condition and qc table."""
    return " ".join(f"--lib {s} {SAMPLES.at[s, 'condition']} {qc}"
                    for s, qc in zip(samples_of(wc.exp, "ribo"), input.qc))


rule rpf_qc_plot:
    input:
        qc=rpf_qc_tables,
        script=workflow.source_path("../scripts/plot_rpf_qc.py"),
    output:
        f"{RPF_QC_OUT}/{{exp}}_rpf_qc.pdf",
    wildcard_constraints:
        exp=EXPERIMENT_RE,
    params:
        lo=config["read_lengths"][0],
        hi=config["read_lengths"][1],
        libraries=rpf_qc_libraries,
    log:
        f"{RPF_QC_LOG}.log",
    conda:
        "../envs/rpf_qc.yaml"
    shell:
        "python {input.script} {output} {wildcards.exp} {params.lo} {params.hi} {params.libraries} > {log} 2>&1"
