"""Blacklist transcripts below an expression floor in the matched total RNA-seq.

Port of autofilter_rnaseq.R (Frederick Korbel, eIF pipeline), same rule: keep
a transcript if TPM >= min_tpm in >= min_samples libraries. Transcripts absent
from a quant.sf count as TPM 0.

Usage: autofilter_rnaseq.py <transcriptome.fa> <blacklist_out> <min_tpm> <min_samples> <quant.sf>...
"""
import sys

import pandas as pd


def main():
    if len(sys.argv) < 6:
        sys.exit(__doc__)
    fasta, blacklist_out, min_tpm, min_samples, *quant_files = sys.argv[1:]
    min_tpm, min_samples = float(min_tpm), int(min_samples)

    with open(fasta) as f:
        txids = [line[1:].split()[0] for line in f if line.startswith(">")]

    tpm = pd.concat(
        [pd.read_csv(q, sep="\t", index_col="Name")["TPM"] for q in quant_files], axis=1
    ).reindex(txids).fillna(0)
    keep = ((tpm >= min_tpm).sum(axis=1) >= min_samples).to_numpy()
    blacklist = [tx for tx, k in zip(txids, keep) if not k]

    with open(blacklist_out, "w") as out:
        out.writelines(tx + "\n" for tx in blacklist)
    print(f"{len(blacklist)} of {len(txids)} transcripts blacklisted "
          f"(TPM >= {min_tpm:g} in >= {min_samples} of {len(quant_files)} libraries required)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
