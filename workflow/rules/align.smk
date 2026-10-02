# Import wf-riboseq-align once, for every experiment in samples.csv (its
# `experiment` column, <target>_<timepoint>): each experiment's outputs land in
# results/align/<experiment>/, and `mode: filtered` builds each one its own
# RNA-seq-filtered reference. What does not depend on the experiment (the
# contaminant index, the unfiltered RNA-seq index) is built once, in
# results/align/.

# wf-riboseq-align at a pinned commit, fetched from GitHub; config
# `align_snakefile` (a local checkout's workflow/Snakefile) replaces it.
ALIGN_SNAKEFILE = config.get("align_snakefile") or github(
    "katarinagresova/wf-riboseq-align", path="workflow/Snakefile", commit="56342d9932b123a49ee2f49f01f21b07352e7ff8"
)


ALIGN_CONFIG = {
    **config["align"],
    # the RNA-seq salmon index's decoys (wf-riboseq-align 7ac127c on); the same genome as ribokit's
    "human_genome_fa": config["references"]["human_genome_fa"],
    "samples": config["samples"],
    "RESULTS_DIR": f"{RESULTS_DIR}/align",
    "LOG_DIR": f"{LOG_DIR}/align",
}


module align:
    snakefile: ALIGN_SNAKEFILE
    config: ALIGN_CONFIG


use rule * from align as align_*
