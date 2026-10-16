#!/usr/bin/env python3
"""Pick the contiguous window of one chromosome with the most back-splice junctions (BSJs).

Only BSJs that lie completely inside a gene that itself lies completely inside the window are
counted, so that the rebased mini GTF still annotates them.

usage: select_region.py <bsj.bed> <genes.gtf> <chrom> <window_size> [step]
prints: <chrom> <start> <end> <n_bsj>     (1-based, inclusive start, inclusive end)
"""
import re
import sys
from collections import defaultdict

bed, gtf, chrom, size = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
step = int(sys.argv[5]) if len(sys.argv) > 5 else 5000

genes = defaultdict(lambda: [10**12, 0])
for line in open(gtf):
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[0] != chrom:
        continue
    gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
    genes[gid][0] = min(genes[gid][0], int(f[3]))
    genes[gid][1] = max(genes[gid][1], int(f[4]))
genes = sorted(genes.values())

bsjs = set()
for line in open(bed):
    f = line.split("\t")
    if f[0] == chrom:
        bsjs.add((int(f[1]) + 1, int(f[2])))  # BED -> 1-based inclusive
bsjs = sorted(bsjs)

chrom_len = max(g[1] for g in genes)
best = (-1, 1)
for start in range(1, chrom_len - size + 2, step):
    end = start + size - 1
    inside = [g for g in genes if g[0] >= start and g[1] <= end]
    n = sum(1 for s, e in bsjs if s >= start and e <= end and any(g[0] <= s and e <= g[1] for g in inside))
    if n > best[0]:
        best = (n, start)
n, start = best
print(chrom, start, start + size - 1, n)
