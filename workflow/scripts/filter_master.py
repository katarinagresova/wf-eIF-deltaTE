"""The filtered master table: per experiment, the transcripts expressed in the reference condition.

Keeps a row of master_table.py's table if, in every library of the reference condition of its experiment, the Ribo-seq
ritpm >= min_ribo_tpm and the RNA-seq TPM >= min_rna_tpm. A missing value fails (a transcript without a CDS has no
ritpm). An experiment's libraries are its <read_type>_tpm_<reference>_rep<rep> columns with a value in any of its rows.
Rows and columns are otherwise the full table's, padj included; the index is the row's number in the full table.

Usage: filter_master.py <master.csv> <out.csv> <reference> <min_ribo_tpm> <min_rna_tpm>
"""
import re
import sys

import pandas as pd


def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    master_csv, out, reference, min_ribo_tpm, min_rna_tpm = sys.argv[1:]
    min_tpm = {"ribo": float(min_ribo_tpm), "rna": float(min_rna_tpm)}

    master = pd.read_csv(master_csv, index_col=0, float_precision="round_trip")   # rows written back unchanged
    keep = pd.Series(True, index=master.index)
    for exp, rows in master.groupby("factor_time", sort=False):
        for read_type, floor in min_tpm.items():
            libraries = [c for c in rows if re.fullmatch(rf"{read_type}_tpm_{re.escape(reference)}_rep\d+", c)
                         and rows[c].notna().any()]
            if not libraries:
                sys.exit(f"{exp}: no {read_type}_tpm_{reference}_rep<rep> column")
            keep.loc[rows.index] &= rows[libraries].ge(floor).all(axis=1)
        print(f"{exp}: {keep.loc[rows.index].sum()} of {len(rows)} rows kept")

    master[keep].to_csv(out)
    print(f"{keep.sum()} of {len(master)} rows -> {out}")


if __name__ == "__main__":
    main()
