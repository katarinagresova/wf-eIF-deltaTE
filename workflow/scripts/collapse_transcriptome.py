"""Collapse a transcriptome by shared CDS: Ribo-seq reads over the CDS cannot
tell apart transcripts whose CDS exon chain is identical, so each such group
(regardless of gene_id -- two gene_ids can be called at the same locus) keeps
one representative: prefer gene_biotype protein_coding, then highest TPM (an
external short-read reference), then the smallest transcript_id. Transcripts
with no CDS at all are grouped by gene_id alone (one representative per gene,
same tie-break) -- they carry no Ribo-seq ambiguity, only the RNA-seq side
needs a representative for them.

Before that, every transcript of a GENCODE readthrough gene goes (RPS10-NUDT3:
every GENCODE transcript of the gene tagged readthrough_transcript, or some
tagged and a hyphenated gene_name). Its CDS spans two genes' CDSs, so the
reads of both are shared with it; matched on the unversioned gene_id.

A second pass then does the same within a gene for near-duplicates: a
transcript >=99% covered (by CDS length) by a same-gene sibling's CDS has
essentially no genomic position of its own to contribute a distinguishing
Ribo-seq read, so clusters connected by that >=99% relation (either
direction) collapse to one representative, same tie-break.

A third pass then keeps one transcript per gene_id: any gene still holding
more than one survivor collapses to one, same tie-break (Frederick's
one-per-gene rule, applied after the CDS passes rather than instead of them,
so a CDS shared across gene_ids is still resolved genome-wide first).

A fourth pass then merges CDSs that share most of their sequence, across
genes and loci (paralogs: SMN1/SMN2, EIF3C/EIF3CL): a CDS of which fewer than
<min_unique_fraction> of the k-mer positions are in no other CDS joins the CDS
it shares most k-mers with, and each connected cluster keeps one, same
tie-break. Repeated until no CDS is below the fraction.

Usage: collapse_transcriptome.py <raw.gtf> <raw.fa> <tpm.csv> <gencode.gtf> <k> <min_unique_fraction>
       <out.gtf> <out.fa> <out_report.tsv>
"""
import csv
import re
import sys
from collections import Counter, defaultdict

import numpy as np

ATTR_RE = re.compile(r'(\w+) "([^"]*)"')


def parse_attrs(field):
    return dict(ATTR_RE.findall(field))


def read_gtf(gtf_path):
    """transcript_id -> metadata (from `transcript` rows), -> CDS and -> exon interval lists."""
    tx_meta = {}
    tx_cds = defaultdict(list)
    tx_exons = defaultdict(list)
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
            elif feature in ("CDS", "exon"):
                tid = parse_attrs(attrs)["transcript_id"]
                (tx_cds if feature == "CDS" else tx_exons)[tid].append((int(start), int(end)))
    return tx_meta, tx_cds, tx_exons


def read_tpm(tpm_path):
    tpm = {}
    with open(tpm_path) as f:
        for row in csv.DictReader(f):
            try:
                tpm[row["transcript_id"]] = float(row["meanTPM"])
            except (KeyError, ValueError):
                pass
    return tpm


def read_fasta(fa_path):
    seqs, name = {}, None
    with open(fa_path) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            else:
                seqs[name].append(line.strip())
    return {t: "".join(parts) for t, parts in seqs.items()}


def read_readthrough_genes(gencode_gtf):
    """Unversioned gene_ids of GENCODE's readthrough genes: every transcript tagged
    readthrough_transcript, or some and a hyphenated gene_name (TIMM23B-AGAP6). A
    gene with one readthrough isoform among normal ones (SGPL1, CARMIL2) is not one."""
    n_tx, n_readthrough, name = Counter(), Counter(), {}
    with open(gencode_gtf) as f:
        for line in f:
            if "\ttranscript\t" not in line:
                continue
            fields = line.split("\t", 8)
            if fields[2] != "transcript":
                continue
            a = parse_attrs(fields[8])  # `tag` repeats, so look for it in the raw field
            gene = a["gene_id"].split(".")[0]
            n_tx[gene] += 1
            n_readthrough[gene] += 'tag "readthrough_transcript"' in fields[8]
            name[gene] = a["gene_name"]
    return {g for g, n in n_readthrough.items() if n and (n == n_tx[g] or "-" in name[g])}


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


CODE = np.full(256, 4, np.uint8)
CODE[np.frombuffer(b"ACGT", np.uint8)] = np.arange(4)


def cds_sequence(t, tx_meta, tx_cds, tx_exons, seq):
    """The CDS of t, 5' to 3' as in the transcript sequence, from its exon and CDS intervals."""
    strand = tx_meta[t]["strand"]
    exons = sorted(tx_exons[t], reverse=strand == "-")
    assert len(seq) == cds_len(exons), f"{t}: sequence length != exon length"
    parts, offset = [], 0
    for s, e in exons:
        for cs, ce in tx_cds[t]:
            if s <= cs and ce <= e:
                a = offset + (cs - s if strand == "+" else e - ce)
                parts.append((a, seq[a:a + ce - cs + 1]))
        offset += e - s + 1
    cds = "".join(p for _, p in sorted(parts))
    assert len(cds) == cds_len(tx_cds[t]), f"{t}: a CDS interval is not within one exon"
    return cds.upper()


def kmers_of(seq, k):
    """The k-mers of seq (A/C/G/T only) as 2-bit integers, one per position."""
    x = CODE[np.frombuffer(seq.encode(), np.uint8)]
    if len(x) < k:
        return np.empty(0, np.uint64)
    w = np.lib.stride_tricks.sliding_window_view(x, k)
    w = w[(w < 4).all(axis=1)].astype(np.uint64)
    return (w << (2 * np.arange(k - 1, -1, -1, dtype=np.uint64))).sum(axis=1)


def low_unique_cds(txs, cds, k, min_unique):
    """The CDSs of txs with < min_unique of their k-mer positions in no other CDS of
    txs -> (fraction unique, the CDS sharing most of its distinct k-mers)."""
    kmers = [kmers_of(cds[t], k) for t in txs]
    ends = np.cumsum([len(x) for x in kmers])
    owner = np.repeat(np.arange(len(txs)), [len(x) for x in kmers])
    kmers = np.concatenate(kmers)
    order = np.lexsort((owner, kmers))
    km, ow = kmers[order], owner[order]
    distinct = np.r_[True, (km[1:] != km[:-1]) | (ow[1:] != ow[:-1])]
    km, ow = km[distinct], ow[distinct]  # each (k-mer, CDS) once, sorted by k-mer
    starts = np.flatnonzero(np.r_[True, km[1:] != km[:-1]])
    n_cds = np.diff(np.r_[starts, len(km)])
    uniq = km[starts]
    n_pos = np.bincount(owner, minlength=len(txs))
    n_unique = np.bincount(owner, weights=n_cds[np.searchsorted(uniq, kmers)] == 1, minlength=len(txs))

    low = {}
    for i in np.flatnonzero((n_pos > 0) & (n_unique < min_unique * n_pos)):
        j = np.searchsorted(uniq, np.unique(kmers[ends[i] - n_pos[i]:ends[i]]))
        j = j[n_cds[j] > 1]
        others = ow[np.concatenate([np.arange(starts[x], starts[x] + n_cds[x]) for x in j])]
        shared = np.bincount(others[others != i], minlength=len(txs))
        low[txs[i]] = (n_unique[i] / n_pos[i], txs[int(shared.argmax())])
    return low


def merge_shared_kmers(whitelist, tx_meta, tx_cds, tx_exons, seqs, tpm, k, min_unique):
    """Across genes, cluster each CDS below min_unique unique k-mer positions with
    the CDS it shares most k-mers with and keep one per cluster; repeated until
    none is below. Returns the surviving set and the (members, winner) of each
    cluster (for the report, group_type "kmer_shared")."""
    cds = {t: cds_sequence(t, tx_meta, tx_cds, tx_exons, seqs[t]) for t in whitelist if t in tx_cds}
    merged = []
    while True:
        txs = sorted(cds)
        low = low_unique_cds(txs, cds, k, min_unique)
        if not low:
            break
        parent = {t: t for t in txs}
        for t, (_, partner) in low.items():
            _union(parent, t, partner)
        clusters = defaultdict(list)
        for t in txs:
            clusters[_find(parent, t)].append(t)
        for members in clusters.values():
            if len(members) < 2:
                continue
            winner = resolve_group(members, tx_meta, tpm)
            merged.append((members, winner))
            for m in members:
                if m != winner:
                    del cds[m]
    dropped = {m for members, winner in merged for m in members if m != winner}
    return whitelist - dropped, merged


def write_report(report_path, groups, tx_meta, tpm):
    """groups: (group_type, members, winner); winner None = every member dropped."""
    header = ["group_type", "chrom", "strand", "n_members", "cross_gene",
              "winner_transcript_id", "winner_gene_id", "winner_gene_name",
              "winner_gene_biotype", "winner_tpm",
              "dropped_transcript_ids", "dropped_gene_ids", "dropped_gene_names",
              "dropped_gene_biotypes", "dropped_tpms"]
    with open(report_path, "w") as f:
        f.write("\t".join(header) + "\n")
        for group_type, members, winner in groups:
            if len(members) < 2 and winner is not None:
                continue
            dropped = [m for m in members if m != winner]
            cross_gene = len(set(tx_meta[m]["gene_id"] for m in members)) > 1
            w = tx_meta[winner] if winner else {}
            locus = tx_meta[winner or members[0]]
            row = [
                group_type, locus["chrom"], locus["strand"], str(len(members)),
                str(cross_gene),
                winner or "", w.get("gene_id", ""), w.get("gene_name", ""),
                w.get("gene_biotype", ""), repr(tpm.get(winner, 0.0)) if winner else "",
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
    if len(sys.argv) != 10:
        sys.exit(__doc__)
    raw_gtf, raw_fa, tpm_path, gencode_gtf, k, min_unique, out_gtf, out_fa, report_path = sys.argv[1:]
    k, min_unique = int(k), float(min_unique)

    tx_meta, tx_cds, tx_exons = read_gtf(raw_gtf)
    seqs = read_fasta(raw_fa)
    assert len(tx_meta) == len(seqs), "transcript rows in the GTF and sequences in the fasta disagree"
    tpm = read_tpm(tpm_path)

    readthrough = read_readthrough_genes(gencode_gtf)
    readthrough_txs = defaultdict(list)
    for t, m in tx_meta.items():
        if m["gene_id"].split(".")[0] in readthrough:
            readthrough_txs[m["gene_id"]].append(t)
    report = [("readthrough", sorted(members), None) for members in readthrough_txs.values()]
    kept = {t: m for t, m in tx_meta.items() if m["gene_id"] not in readthrough_txs}

    groups = build_groups(kept, tx_cds)
    stage1_winners = [resolve_group(members, tx_meta, tpm) for _, members in groups]
    report += [(group_type, members, w) for (group_type, members), w in zip(groups, stage1_winners)]
    cds_winners = [w for (group_type, _), w in zip(groups, stage1_winners) if group_type == "cds"]
    nocds_winners = [w for (group_type, _), w in zip(groups, stage1_winners) if group_type == "nocds"]

    overlap_clusters = merge_near_duplicates(cds_winners, tx_meta, tx_cds)
    overlap_winners = [resolve_group(members, tx_meta, tpm) for members in overlap_clusters]
    report += [("cds_overlap99", members, w) for members, w in zip(overlap_clusters, overlap_winners)]

    survivors = set(overlap_winners) | set(nocds_winners)
    gene_clusters, whitelist = force_one_per_gene(survivors, tx_meta, tpm)
    report += [("gene_forced", members, resolve_group(members, tx_meta, tpm)) for members in gene_clusters]

    whitelist, merged = merge_shared_kmers(whitelist, tx_meta, tx_cds, tx_exons, seqs, tpm, k, min_unique)
    report += [("kmer_shared", members, w) for members, w in merged]

    write_report(report_path, report, tx_meta, tpm)
    write_gtf(raw_gtf, out_gtf, whitelist)
    write_fasta(raw_fa, out_fa, whitelist)


if __name__ == "__main__":
    main()
