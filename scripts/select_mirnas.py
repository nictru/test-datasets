#!/usr/bin/env python3
"""Pick miRNAs with the most circRNA binding sites and write mature_mini.fa + a synthetic expression table.

Ranking: number of circular transcripts (ids like chrI:start-end:strand) predicted by BOTH miranda and
targetscan (so that majority_vote works with min_tools 1 and 2); ties are broken by name.

usage: select_mirnas.py <binding_sites_dir> <mature.fa> <n> <out_mature.fa> <out_expression.tsv> <sample> [<sample> ...]
"""
import glob
import random
import sys
from collections import defaultdict

sites_dir, mature, n, out_fa, out_tsv, *samples = sys.argv[1:]
n = int(n)

def pairs(tool):
    d = defaultdict(set)
    for path in glob.glob(f"{sites_dir}/tools/{tool}/unified/*.tsv"):
        for line in open(path):
            f = line.split("\t")
            if len(f) > 1 and f[1].startswith("chr"):
                d[f[0]].add(f[1])
    return d

mi, ts = pairs("miranda"), pairs("targetscan")
score = {m: len(mi[m] & ts[m]) for m in mi if m in ts}
chosen = sorted(score, key=lambda m: (-score[m], m))[:n]

seqs, name = {}, None
for line in open(mature):
    if line.startswith(">"):
        name = line[1:].split()[0]
        seqs[name] = [line.rstrip("\n"), ""]
    else:
        seqs[name][1] += line.strip()
with open(out_fa, "w") as o:
    for m in sorted(chosen):
        o.write(seqs[m][0] + "\n" + seqs[m][1] + "\n")

rng = random.Random(42)
with open(out_tsv, "w") as o:
    o.write("miRNA\t" + "\t".join(samples) + "\n")
    for m in sorted(chosen):
        base = int(rng.lognormvariate(5, 1)) + 20
        o.write(m + "\t" + "\t".join(str(max(10, int(base * rng.uniform(0.6, 1.6)))) for _ in samples) + "\n")
print({m: score[m] for m in sorted(chosen)}, file=sys.stderr)
