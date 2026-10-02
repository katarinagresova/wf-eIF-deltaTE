"""Collapse a transcriptome by shared CDS: Ribo-seq reads over the CDS cannot
tell apart transcripts whose CDS exon chain is identical, so each such group
(regardless of gene_id -- two gene_ids can be called at the same locus) keeps
one representative: prefer gene_biotype protein_coding, then highest TPM (an
external short-read reference), then the smallest transcript_id. Transcripts
with no CDS at all are grouped by gene_id alone (one representative per gene,
same tie-break) -- they carry no Ribo-seq ambiguity, only the RNA-seq side
needs a representative for them.

A second pass then does the same within a gene for near-duplicates: a
transcript >=99% covered (by CDS length) by a same-gene sibling's CDS has
essentially no genomic position of its own to contribute a distinguishing
Ribo-seq read, so clusters connected by that >=99% relation (either
direction) collapse to one representative, same tie-break.

A third pass then keeps one transcript per gene_id: any gene still holding
more than one survivor collapses to one, same tie-break (Frederick's
one-per-gene rule, applied after the CDS passes rather than instead of them,
so a CDS shared across gene_ids is still resolved genome-wide first).

Usage: collapse_transcriptome.py <raw.gtf> <raw.fa> <tpm.csv> <out.gtf> <out.fa> <out_report.tsv>
"""
import csv
import re
import sys
from collections import defaultdict

ATTR_RE = re.compile(r'(\w+) "([^"]*)"')


def parse_attrs(field):
    return dict(ATTR_RE.findall(field))


def read_gtf(gtf_path):
    """transcript_id -> metadata (from `transcript` rows) and -> CDS interval list."""
    tx_meta = {}
    tx_cds = defaultdict(list)
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = line.rstrip("\n").split("\t")
            if feature == "transcript":
                a = parse_attrs(attrs)
                tx_meta[a["transcript_id"]] = {
                    "chrom": chrom,
                    "strand": strand,
                    "gene_id": a["gene_id"],
                    "gene_name": a["gene_name"],
                    "gene_biotype": a["gene_biotype"],
                    "transcript_biotype": a["transcript_biotype"],
                }
            elif feature == "CDS":
                tid = parse_attrs(attrs)["transcript_id"]
                tx_cds[tid].append((int(start), int(end)))
    return tx_meta, tx_cds


def read_tpm(tpm_path):
    tpm = {}
    with open(tpm_path) as f:
        for row in csv.DictReader(f):
            try:
                tpm[row["transcript_id"]] = float(row["meanTPM"])
            except (KeyError, ValueError):
                pass
    return tpm


def build_groups(tx_meta, tx_cds):
    """list of (group_type, members) -- CDS-bearing transcripts grouped genome-wide
    by (chrom, strand, CDS intervals); the rest grouped by gene_id."""
    cds_groups = defaultdict(list)
    nocds_groups = defaultdict(list)
    for tid, meta in tx_meta.items():
        if tid in tx_cds:
            sig = (meta["chrom"], meta["strand"], tuple(sorted(tx_cds[tid])))
            cds_groups[sig].append(tid)
        else:
            nocds_groups[meta["gene_id"]].append(tid)
    groups = [("cds", members) for members in cds_groups.values()]
    groups += [("nocds", members) for members in nocds_groups.values()]
    return groups


def resolve_group(members, tx_meta, tpm):
    """The representative transcript_id: protein_coding gene_biotype preferred,
    then highest TPM, then smallest transcript_id."""
    protein_coding = [m for m in members if tx_meta[m]["gene_biotype"] == "protein_coding"]
    pool = protein_coding if protein_coding else members
    return min(pool, key=lambda m: (-tpm.get(m, 0.0), m))


OVERLAP_THRESHOLD = 0.99


def cds_len(intervals):
    return sum(e - s + 1 for s, e in intervals)


def overlap_len(a, b):
    """total overlap length between two sorted, disjoint interval lists."""
    i, j, total = 0, 0, 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if lo <= hi:
            total += hi - lo + 1
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total


def _find(parent, x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def _union(parent, a, b):
    ra, rb = _find(parent, a), _find(parent, b)
    if ra != rb:
        parent[ra] = rb


def merge_near_duplicates(cds_winners, tx_meta, tx_cds):
    """Within a gene, cluster transcripts connected by >=OVERLAP_THRESHOLD
    CDS-length coverage (either direction) and keep one representative per
    cluster. Returns the surviving set and the clusters of size > 1 (for the
    report, group_type "cds_overlap99")."""
    by_gene = defaultdict(list)
    for t in cds_winners:
        by_gene[tx_meta[t]["gene_id"]].append(t)

    parent = {t: t for t in cds_winners}
    for txs in by_gene.values():
        for i in range(len(txs)):
            for j in range(i + 1, len(txs)):
                t1, t2 = txs[i], txs[j]
                ov = overlap_len(tx_cds[t1], tx_cds[t2])
                if ov == 0:
                    continue
                if ov / cds_len(tx_cds[t1]) >= OVERLAP_THRESHOLD or ov / cds_len(tx_cds[t2]) >= OVERLAP_THRESHOLD:
                    _union(parent, t1, t2)

    clusters = defaultdict(list)
    for t in cds_winners:
        clusters[_find(parent, t)].append(t)
    return list(clusters.values())


def force_one_per_gene(survivors, tx_meta, tpm):
    """Group the survivors by gene_id and keep one representative per gene.
    Returns the gene clusters (for the report, group_type "gene_forced") and
    the final whitelist."""
    by_gene = defaultdict(list)
    for t in sorted(survivors):
        by_gene[tx_meta[t]["gene_id"]].append(t)
    clusters = list(by_gene.values())
    whitelist = {resolve_group(members, tx_meta, tpm) for members in clusters}
    return clusters, whitelist


def write_report(report_path, groups, tx_meta, tpm):
    header = ["group_type", "chrom", "strand", "n_members", "cross_gene",
              "winner_transcript_id", "winner_gene_id", "winner_gene_name",
              "winner_gene_biotype", "winner_tpm",
              "dropped_transcript_ids", "dropped_gene_ids", "dropped_gene_names",
              "dropped_gene_biotypes", "dropped_tpms"]
    with open(report_path, "w") as f:
        f.write("\t".join(header) + "\n")
        for group_type, members in groups:
            if len(members) < 2:
                continue
            winner = resolve_group(members, tx_meta, tpm)
            dropped = [m for m in members if m != winner]
            cross_gene = len(set(tx_meta[m]["gene_id"] for m in members)) > 1
            row = [
                group_type, tx_meta[winner]["chrom"], tx_meta[winner]["strand"], str(len(members)), str(cross_gene),
                winner, tx_meta[winner]["gene_id"], tx_meta[winner]["gene_name"],
                tx_meta[winner]["gene_biotype"], repr(tpm.get(winner, 0.0)),
                ",".join(dropped),
                ",".join(tx_meta[m]["gene_id"] for m in dropped),
                ",".join(tx_meta[m]["gene_name"] for m in dropped),
                ",".join(tx_meta[m]["gene_biotype"] for m in dropped),
                ",".join(repr(tpm.get(m, 0.0)) for m in dropped),
            ]
            f.write("\t".join(row) + "\n")


def write_gtf(raw_gtf, out_gtf, whitelist):
    with open(raw_gtf) as fin, open(out_gtf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            attrs = line.rstrip("\n").split("\t")[8]
            if parse_attrs(attrs)["transcript_id"] in whitelist:
                fout.write(line)


def write_fasta(raw_fa, out_fa, whitelist):
    with open(raw_fa) as fin, open(out_fa, "w") as fout:
        keep = False
        for line in fin:
            if line.startswith(">"):
                keep = line[1:].split()[0] in whitelist
            if keep:
                fout.write(line)


def main():
    if len(sys.argv) != 7:
        sys.exit(__doc__)
    raw_gtf, raw_fa, tpm_path, out_gtf, out_fa, report_path = sys.argv[1:]

    tx_meta, tx_cds = read_gtf(raw_gtf)
    assert len(tx_meta) == sum(1 for line in open(raw_fa) if line.startswith(">")), \
        "transcript rows in the GTF and sequences in the fasta disagree"
    tpm = read_tpm(tpm_path)
    groups = build_groups(tx_meta, tx_cds)

    stage1_winners = [resolve_group(members, tx_meta, tpm) for _, members in groups]
    cds_winners = [w for (group_type, _), w in zip(groups, stage1_winners) if group_type == "cds"]
    nocds_winners = [w for (group_type, _), w in zip(groups, stage1_winners) if group_type == "nocds"]

    overlap_clusters = merge_near_duplicates(cds_winners, tx_meta, tx_cds)
    groups += [("cds_overlap99", members) for members in overlap_clusters if len(members) > 1]
    final_cds_winners = {resolve_group(members, tx_meta, tpm) for members in overlap_clusters}

    survivors = final_cds_winners | set(nocds_winners)
    gene_clusters, whitelist = force_one_per_gene(survivors, tx_meta, tpm)
    groups += [("gene_forced", members) for members in gene_clusters if len(members) > 1]

    write_report(report_path, groups, tx_meta, tpm)
    write_gtf(raw_gtf, out_gtf, whitelist)
    write_fasta(raw_fa, out_fa, whitelist)


if __name__ == "__main__":
    main()
