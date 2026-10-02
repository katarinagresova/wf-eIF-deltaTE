"""Drop blacklisted transcripts from a salmon quant.sf.

The eIF pipeline's results/salmon/data_filtered/<sample>_quant.sf (Frederick Korbel): rows whose Name is in the
blacklist are removed (exact match), the rest sorted by Name. Values are salmon's, copied as written: TPM is not
renormalised over the kept rows, and the reads salmon gave the dropped transcripts are not redistributed. (His R
rewrote the numbers with readr; numerically identical on all 6 eIF4E 4 h libraries, same 15,508 rows in the same
order.)

Usage: filter_quant.py <quant.sf> <blacklist> <quant_out.sf>
"""
import sys


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    quant, blacklist, out = sys.argv[1:]

    with open(blacklist) as f:
        remove = {line.strip() for line in f if line.strip()}
    with open(quant) as f:
        header, *rows = f.readlines()
    kept = sorted((r for r in rows if r.split("\t", 1)[0] not in remove), key=lambda r: r.split("\t", 1)[0])
    with open(out, "w") as f:
        f.writelines([header] + kept)

    print(f"{len(rows) - len(kept)} of {len(rows)} transcripts blacklisted, {len(kept)} kept in {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
