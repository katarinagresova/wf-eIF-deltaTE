"""Length and GC content of every transcript in a fasta.

Usage: transcript_features.py <transcriptome.fa> <out.tsv>
Writes Name, length, gc (fraction of G + C), one row per transcript, in fasta order.
"""
import sys


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    fasta, out = sys.argv[1:]

    seq, name = {}, None
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seq[name] = []
            else:
                seq[name].append(line.strip().upper())
    with open(out, "w") as f:
        f.write("Name\tlength\tgc\n")
        for name, parts in seq.items():
            s = "".join(parts)
            f.write(f"{name}\t{len(s)}\t{(s.count('G') + s.count('C')) / len(s)!r}\n")


if __name__ == "__main__":
    main()
