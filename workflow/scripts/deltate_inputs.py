"""One experiment's inputs to DESeq2 (workflow/rules/deltate.smk, workflow/scripts/deseq2.R).

The Ribo-seq and RNA-seq libraries as two transcript x library matrices, Ribo-seq libraries first, as tximport
(type "salmon", txOut, countsFromAbundance "no") would give them: counts = NumReads, length = EffectiveLength. The rows
are the RNA-seq quants' transcripts (every RNA-seq quant has the same), sorted by Name; a transcript a Ribo-seq quant
lacks (no CDS) gets 0 reads and length 1, and a missing NumReads counts 0. Per Ribo-seq library its spike-in reads
(deseq2.R's spike-in size factors are each library's over the first's, so the order of the --ribo libraries matters).

Usage: deltate_inputs.py <dir> --ribo <sample> <condition> <human quant> <spike-in quant> ...
                               --rna <sample> <condition> <quant> ...
  human quant, spike-in quant  ribokit's quants (Name, Length, EffectiveLength, ritpm, NumReads)
  quant                        salmon's quant.sf
Writes into <dir>:
  sampleTable.tsv             sampleName, assay (ribo / total), condition
  counts.tsv, length.tsv      Name + one column per library
  spike_in_reads.tsv          sampleName, reads (Ribo-seq libraries only)
"""
import argparse
import math
import os

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("dir")
parser.add_argument("--ribo", nargs=4, action="append", required=True,
                    metavar=("SAMPLE", "CONDITION", "HUMAN_QUANT", "SPIKE_IN_QUANT"))
parser.add_argument("--rna", nargs=3, action="append", required=True, metavar=("SAMPLE", "CONDITION", "QUANT"))


def read_quant(path):
    q = pd.read_csv(path, sep="\t", index_col="Name", float_precision="round_trip")
    assert q.index.is_unique, path
    return q


def main():
    args = parser.parse_args()
    ribo = {s: read_quant(path) for s, _, path, _ in args.ribo}
    rna = {s: read_quant(path) for s, _, path in args.rna}
    assert len(ribo) == len(args.ribo) and len(rna) == len(args.rna) and not set(ribo) & set(rna), "duplicate sample"

    names = sorted(next(iter(rna.values())).index)
    for s, q in rna.items():
        assert set(q.index) == set(names), f"{s}: other transcripts than {args.rna[0][0]}"
    for s, q in ribo.items():
        assert set(q.index) <= set(names), f"{s}: transcripts without RNA-seq: {sorted(set(q.index) - set(names))[:5]}"

    quants = ribo | rna
    counts = pd.DataFrame({s: q["NumReads"].reindex(names).fillna(0) for s, q in quants.items()})
    length = pd.DataFrame({s: q["EffectiveLength"].reindex(names).fillna(1 if s in ribo else math.nan)
                           for s, q in quants.items()}).astype(float)
    assert (length > 0).all().all(), "a length is missing or not positive"

    spike_in = pd.Series({s: math.fsum(read_quant(path)["NumReads"].dropna()) for s, _, _, path in args.ribo})

    pd.DataFrame({"sampleName": list(quants),
                  "assay": ["ribo"] * len(ribo) + ["total"] * len(rna),
                  "condition": [c for _, c, *_ in args.ribo] + [c for _, c, _ in args.rna]}).to_csv(
        os.path.join(args.dir, "sampleTable.tsv"), sep="\t", index=False)
    counts.to_csv(os.path.join(args.dir, "counts.tsv"), sep="\t", index_label="Name")
    length.to_csv(os.path.join(args.dir, "length.tsv"), sep="\t", index_label="Name")
    spike_in.rename_axis("sampleName").rename("reads").to_csv(os.path.join(args.dir, "spike_in_reads.tsv"), sep="\t")

    print(f"{len(names)} transcripts, {len(set().union(*(q.index for q in ribo.values())))} in a Ribo-seq quant; "
          f"{len(ribo)} Ribo-seq + {len(rna)} RNA-seq libraries")
    print(f"spike-in reads: {spike_in.to_dict()}")


if __name__ == "__main__":
    main()
