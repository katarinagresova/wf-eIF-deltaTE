"""Drop blacklisted transcripts from a GTF.

The GTF half of wf-riboseq-align's filter_transcriptome.py (port of
make_resources_autofilter.sh, Frederick Korbel, eIF pipeline). Matches
transcript ids EXACTLY: the original's `grep -F` is a substring match, which
over-matches versioned ids (ENST...123.1 is inside ENST...123.10). GTF lines
without a transcript_id (gene lines, comments) are kept.

Usage: filter_gtf.py <gtf_in> <blacklist> <gtf_out>
"""
import re
import sys

TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    gtf_in, blacklist, gtf_out = sys.argv[1:]

    with open(blacklist) as f:
        remove = {line.strip() for line in f if line.strip()}

    kept = set()
    with open(gtf_in) as src, open(gtf_out, "w") as out:
        for line in src:
            match = TRANSCRIPT_ID.search(line)
            if match is None or match.group(1) not in remove:
                out.write(line)
                if match is not None:
                    kept.add(match.group(1))

    print(f"{len(remove)} transcripts blacklisted, {len(kept)} kept in {gtf_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
