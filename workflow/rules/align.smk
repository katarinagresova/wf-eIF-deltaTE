# Import wf-riboseq-align once, for every experiment in samples.csv (its
# `experiment` column, <target>_<timepoint>): each experiment's outputs land in
# results/align/<experiment>/, and `mode: filtered` builds each one its own
# RNA-seq-filtered reference. What does not depend on the experiment (the
# contaminant index, the unfiltered RNA-seq index) is built once, in
# results/align/.

# wf-riboseq-align at a pinned commit, fetched from GitHub; config
# `align_snakefile` (a local checkout's workflow/Snakefile) replaces it.
ALIGN_SNAKEFILE = config.get("align_snakefile") or github(
    "katarinagresova/wf-riboseq-align", path="workflow/Snakefile", commit="7d108a52ec2285821ca5008e50186f55b093774e"
)


ALIGN_CONFIG = {
    **config["align"],
    "samples": config["samples"],
    "RESULTS_DIR": f"{RESULTS_DIR}/align",
    "LOG_DIR": f"{LOG_DIR}/align",
}


module align:
    snakefile: ALIGN_SNAKEFILE
    config: ALIGN_CONFIG


use rule * from align as align_*
