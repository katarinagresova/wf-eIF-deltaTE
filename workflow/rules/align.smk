# Import wf-riboseq-align once, for every experiment in samples.csv (its
# `experiment` column, <target>_<timepoint>): it aligns the Ribo-seq only, each
# experiment's outputs in results/align/<experiment>/. What does not depend on
# the experiment (the Ribo-seq reference, the contaminant and transcriptome STAR
# indexes) is built once, in results/align/. The RNA-seq: rnaseq.smk.

from snakemake.io import sourcecache_entry

# wf-riboseq-align at a pinned commit, fetched from GitHub; config
# `align_snakefile` (a local checkout's workflow/Snakefile) replaces it.
ALIGN_REPO, ALIGN_COMMIT = "katarinagresova/wf-riboseq-align", "7c5cee5c9d22528135342a394fc03a7ef86d28d5"
ALIGN_SNAKEFILE = config.get("align_snakefile") or github(ALIGN_REPO, path="workflow/Snakefile", commit=ALIGN_COMMIT)

# The Ribo-seq contaminant set wf-riboseq-align ships (its build_contaminants
# output), from the same commit, through snakemake's source cache as the
# module's own scripts are; config `align: contaminants_fa` replaces it.
_CONTAMINANTS = github(ALIGN_REPO, path="resources/contaminants_built.fa", commit=ALIGN_COMMIT)


ALIGN_CONFIG = {
    "contaminants_fa": sourcecache_entry(workflow.sourcecache.get_path(_CONTAMINANTS),
                                         _CONTAMINANTS.get_path_or_uri()),
    **config["align"],
    # our own CDS-based collapse (collapse_transcriptome.smk), not config's
    "human_transcriptome_fa": COLLAPSED_FA,
    "samples": config["samples"],
    "RESULTS_DIR": f"{RESULTS_DIR}/align",
    "LOG_DIR": f"{LOG_DIR}/align",
}


module align:
    snakefile: ALIGN_SNAKEFILE
    config: ALIGN_CONFIG


use rule * from align as align_*
