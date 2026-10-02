# Re-collapse the raw HCT116 transcriptome ourselves (Frederick's own collapsed
# reference picks one representative per gene_id, biotype-blind, and never
# checks for a CDS shared ACROSS gene_ids): a group of transcripts with an
# identical CDS exon chain (same chrom/strand/CDS intervals) is indistinguishable
# to a Ribo-seq read over the CDS, so each group keeps exactly one
# representative (scripts/collapse_transcriptome.py says the tie-break).
# Transcripts with no CDS at all are grouped by gene_id alone. Then same-gene
# CDS near-duplicates (>=99% covered) merge, and last every gene keeps one
# transcript (validation/one_per_gene/: the `hybrid` arm). See
# results/resources/collapse_report.tsv for every group resolved.

_COLLAPSE_CONFIG = config["collapse_transcriptome"]
RAW_TRANSCRIPTOME_FA = _COLLAPSE_CONFIG["raw_fa"]
RAW_TRANSCRIPTOME_GTF = _COLLAPSE_CONFIG["raw_gtf"]
TPM_TABLE = _COLLAPSE_CONFIG["tpm_table"]

COLLAPSE_DIR = f"{RESULTS_DIR}/resources"
COLLAPSED_FA = f"{COLLAPSE_DIR}/HCT116_Txome_WT.v1.1_nonAUG.sort.collapsed.fa"
COLLAPSED_GTF = f"{COLLAPSE_DIR}/HCT116_Txome_WT.v1.1_nonAUG.sort.collapsed.gtf"
COLLAPSE_REPORT = f"{COLLAPSE_DIR}/collapse_report.tsv"


rule collapse_transcriptome:
    input:
        gtf=RAW_TRANSCRIPTOME_GTF,
        fa=RAW_TRANSCRIPTOME_FA,
        tpm=TPM_TABLE,
        script=workflow.source_path("../scripts/collapse_transcriptome.py"),
    output:
        gtf=COLLAPSED_GTF,
        fa=COLLAPSED_FA,
        report=COLLAPSE_REPORT,
    log:
        f"{LOG_DIR}/collapse_transcriptome.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.gtf} {input.fa} {input.tpm} {output.gtf} {output.fa} {output.report} 2> {log}"
