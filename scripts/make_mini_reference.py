#!/usr/bin/env python3
"""Cut a region out of a genome FASTA/GTF and rebase it to start at 1 (chromosome name is kept).

Only genes that lie completely inside the region are kept.

usage: make_mini_reference.py <fa> <gtf> <chrom> <start> <end> <out.fa> <out.gtf>
"""
import re
import sys
from collections import defaultdict

fa, gtf, chrom, start, end, out_fa, out_gtf = sys.argv[1:8]
start, end = int(start), int(end)

seq, take = [], False
for line in open(fa):
    if line.startswith(">"):
        take = line[1:].split()[0] == chrom
    elif take:
        seq.append(line.strip())
seq = "".join(seq)[start - 1:end]
with open(out_fa, "w") as o:
    o.write(f">{chrom}\n")
    for i in range(0, len(seq), 60):
        o.write(seq[i:i + 60] + "\n")

rows, span = [], defaultdict(lambda: [10**12, 0])
for line in open(gtf):
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[0] != chrom:
        continue
    gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
    span[gid][0] = min(span[gid][0], int(f[3]))
    span[gid][1] = max(span[gid][1], int(f[4]))
    rows.append((gid, f))
keep = {g for g, (s, e) in span.items() if s >= start and e <= end}
n = 0
with open(out_gtf, "w") as o:
    for gid, f in rows:
        if gid in keep:
            f[3], f[4] = str(int(f[3]) - start + 1), str(int(f[4]) - start + 1)
            o.write("\t".join(f) + "\n")
            n += 1
print(f"{out_fa}: {len(seq)} bp; {out_gtf}: {len(keep)} genes, {n} features", file=sys.stderr)
