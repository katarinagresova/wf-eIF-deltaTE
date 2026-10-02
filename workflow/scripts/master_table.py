"""The eIF master table (the paper's TE table, read by wf-transmod-data).

A port of Frederick Korbel's notebooks/compile_ribopipe_results_autofilter.ipynb, the cells that wrote
results/20260324-eif-master.csv. Logic and columns unchanged:
- per experiment, an outer join on transcript of log2FC/padj from deseq_res_deltaTE_yeastnorm (TE),
  deseq_res_diffribo_yeastnorm (RPF) and deseq_res_difftotal_autonorm (RNA), the ritpm per Ribo-seq library
  and salmon's TPM per RNA-seq library (after the RNA-seq library correction), as <read_type>_tpm_<condition>_rep<rep>;
- the experiments stacked, with factor_time = the experiment name in lower case (eif4e_4h);
- transcript annotation from the GTF's transcript lines (gene_id only where GENCODE_transcript_id is set, as his);
- log2TE_<condition>_rep<rep> = log2(ribo TPM / RNA TPM) where the ribo TPM exists and the RNA TPM is > 0.
The notebook's other cells only print or plot.

Usage: master_table.py <gtf> <out.csv> --experiment <name> <TE.tsv> <RPF.tsv> <RNA.tsv> ...
                                       --library <experiment> <read_type> <condition> <rep> <quant> ...
--experiment: per experiment its deseq_res_{deltaTE_yeastnorm,diffribo_yeastnorm,difftotal_autonorm}.tsv, written
in the order given. --library: per library samples.csv's read_type (ribo: ribokit's human quant; rna: salmon's
quant.sf), condition and rep. The library columns are sorted by condition, read_type, rep, which gives his order
(ribo minusAux, RNA minusAux, ribo plusAux, RNA plusAux); log2TE pairs ribo and RNA by condition and rep.
"""
import argparse
import re
import sys

import numpy as np
import pandas as pd

TPM_COLUMN = {"ribo": "ritpm", "rna": "TPM"}

parser = argparse.ArgumentParser()
parser.add_argument("gtf")
parser.add_argument("out")
parser.add_argument("--experiment", nargs=4, action="append", required=True,
                    metavar=("NAME", "TE", "RPF", "RNA"))
parser.add_argument("--library", nargs=5, action="append", required=True,
                    metavar=("EXPERIMENT", "READ_TYPE", "CONDITION", "REP", "QUANT"))
args = parser.parse_args()
gtf, out = args.gtf, args.out


def library_key(lib):
    _, read_type, condition, rep, _ = lib
    return condition, read_type, int(rep)


dataframes = []
for exp, *tables in args.experiment:
    dfs_to_merge = [pd.read_csv(path, sep="\t", usecols=["Name", "log2FoldChange", "padj"]).rename(
        columns={"log2FoldChange": f"log2FC_{x}", "padj": f"padj_{x}"}) for x, path in zip(("TE", "RPF", "RNA"), tables)]
    libraries = sorted((lib for lib in args.library if lib[0] == exp), key=library_key)
    columns = [f"{read_type}_tpm_{condition}_rep{rep}" for _, read_type, condition, rep, _ in libraries]
    if len(set(columns)) < len(columns):
        sys.exit(f"{exp}: two libraries share read_type, condition and rep: {sorted(columns)}")
    dfs_to_merge += [pd.read_csv(quant, sep="\t", usecols=["Name", TPM_COLUMN[read_type]]).rename(
        columns={TPM_COLUMN[read_type]: column}) for (_, read_type, _, _, quant), column in zip(libraries, columns)]
    merged = dfs_to_merge[0]
    for df in dfs_to_merge[1:]:
        merged = merged.merge(df, on="Name", how="outer")
    merged["factor_time"] = exp.lower()
    dataframes.append(merged)
stacked_df = pd.concat(dataframes, axis=0, ignore_index=True)
if stacked_df["Name"].isna().any():
    sys.exit("transcripts without a Name")

txome = pd.read_csv(gtf, sep="\t", header=None, comment="#",
                    names=["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"])


def parse_attributes(colstring):
    transcript_id = re.search(r'transcript_id "([^"]+)"', colstring)
    gencode_transcript_id = re.search(r'GENCODE_transcript_id "([^"]+)"', colstring)
    gene_id = re.search(r'gene_id "([^"]+)"', colstring)
    structural_category = re.search(r'structural_category "([^"]+)"', colstring)
    transcript_biotype = re.search(r'transcript_biotype "([^"]+)"', colstring)
    transcript_name = re.search(r'transcript_name "([^"]+)"', colstring)
    gene_name = re.search(r'gene_name "([^"]+)"', colstring)
    return {
        "transcript_id": transcript_id.group(1) if transcript_id else "",
        "GENCODE_transcript_id": gencode_transcript_id.group(1) if gencode_transcript_id else "",
        "gene_id": gene_id.group(1) if gencode_transcript_id else "",   # sic, as in the notebook
        "gene_name": gene_name.group(1) if gene_name else "",
        "structural_category": structural_category.group(1) if structural_category else "",
        "transcript_biotype": transcript_biotype.group(1) if transcript_biotype else "",
        "transcript_name": transcript_name.group(1) if transcript_name else "",
    }


transcripts = txome[txome["feature"] == "transcript"].reset_index()["attributes"].apply(parse_attributes)
annot = pd.DataFrame(transcripts.to_list()).rename(columns={"transcript_id": "Name"})
master = stacked_df.merge(right=annot, on="Name", how="left")

# his column order: rep, then condition
pairs = [{(condition, rep) for _, t, condition, rep, _ in args.library if t == read_type} for read_type in ("ribo", "rna")]
with np.errstate(divide="ignore", invalid="ignore"):
    for condition, rep in sorted(pairs[0] & pairs[1], key=lambda cr: (int(cr[1]), cr[0])):
        ribo_col, rna_col = f"ribo_tpm_{condition}_rep{rep}", f"rna_tpm_{condition}_rep{rep}"
        valid = master[ribo_col].notna() & (master[rna_col] > 0)
        master[f"log2TE_{condition}_rep{rep}"] = np.where(valid, np.log2(master[ribo_col] / master[rna_col]), np.nan)

master.to_csv(out)
print(f"{len(master)} rows, {master.factor_time.nunique()} experiments -> {out}")
