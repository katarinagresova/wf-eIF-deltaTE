"""Drop blacklisted transcripts from a transcriptome BAM.

Removes their @SQ lines and every alignment to them (exact name match, as filter_gtf.py). The other
references keep their order, so a coordinate-sorted BAM stays sorted. Records are copied as they
are: FLAG, NH and HI are not updated, so a read can be left with secondary alignments only, and NH
still counts the removed ones. Single-end reads only. Writes the index (<bam_out>.bai) too.

Usage: filter_bam.py <bam_in> <blacklist> <bam_out>
"""
import sys

import pysam


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    bam_in, blacklist, bam_out = sys.argv[1:]

    with open(blacklist) as f:
        remove = {line.strip() for line in f if line.strip()}

    n_aln = n_dropped = n_primary_dropped = 0
    with pysam.AlignmentFile(bam_in, "rb") as src:
        header = src.header.to_dict()
        header["SQ"] = [sq for sq in header["SQ"] if sq["SN"] not in remove]
        new_id = {sq["SN"]: i for i, sq in enumerate(header["SQ"])}
        new_tid = [new_id.get(r, -1) for r in src.references]
        n_found = len(remove & set(src.references))
        with pysam.AlignmentFile(bam_out, "wb", header=header) as out:
            for a in src.fetch(until_eof=True):
                if a.is_paired:
                    sys.exit(f"{bam_in}: {a.query_name} is paired; single-end reads only")
                n_aln += 1
                if a.reference_id >= 0:
                    if new_tid[a.reference_id] < 0:
                        n_dropped += 1
                        n_primary_dropped += not (a.is_secondary or a.is_supplementary)
                        continue
                    a.reference_id = new_tid[a.reference_id]
                out.write(a)
    pysam.index(bam_out)

    print(f"{n_found} of {len(remove)} blacklisted transcripts in the header, removed; "
          f"{n_dropped} of {n_aln} alignments dropped, {n_primary_dropped} of them primary; {bam_out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
