#!/usr/bin/env python3
"""Select read pairs for the mini fastqs.

A pair is kept if (a) both mates are aligned inside the region (from the unsorted STAR BAM, via `samtools view`)
or (b) STAR reported it as a chimeric (BSJ) read with both breakpoints inside the region.
All chimeric (BSJ) pairs are always kept; the remaining budget up to <max_pairs> is filled with the
aligned pairs that have the smallest md5(read name) (deterministic pseudo-random subsample).

usage: pick_read_pairs.py <bam> <chimeric.junction> <chrom> <start> <end> <fq1> <fq2> <out1> <out2> <max_pairs>
"""
import gzip
import hashlib
import io
import subprocess
import sys

bam, junc, chrom, start, end, fq1, fq2, out1, out2, max_pairs = sys.argv[1:11]
start, end, max_pairs = int(start), int(end), int(max_pairs)

def inside(p):
    return start <= p <= end

names = set()
view = subprocess.run(["samtools", "view", "-F", "0x904", bam], capture_output=True, text=True, check=True)
for line in view.stdout.splitlines():
    f = line.split("\t")
    if f[2] == chrom and f[6] in ("=", chrom) and inside(int(f[3])) and inside(int(f[7])):
        names.add(f[0])
aligned = set(names)
chimeric = set()
for line in open(junc):
    f = line.split("\t")
    if len(f) > 10 and f[0] == chrom and f[3] == chrom and inside(int(f[1])) and inside(int(f[4])):
        chimeric.add(f[9])
rest = sorted(aligned - chimeric, key=lambda n: hashlib.md5(n.encode()).hexdigest())
keep = chimeric | set(rest[:max(0, max_pairs - len(chimeric))])
print(f"{len(aligned)} aligned pairs in region, {len(chimeric)} chimeric reads, keeping {len(keep)}", file=sys.stderr)

def records(path):
    with gzip.open(path, "rt") as fh:
        while True:
            rec = [fh.readline() for _ in range(4)]
            if not rec[0]:
                return
            yield rec

n = 0
def open_gz(path):  # mtime=0 keeps the output byte-identical between runs
    return io.TextIOWrapper(gzip.GzipFile(path, "wb", compresslevel=9, mtime=0))

with open_gz(out1) as o1, open_gz(out2) as o2:
    for r1, r2 in zip(records(fq1), records(fq2)):
        if r1[0][1:].split()[0].split("/")[0] in keep:
            o1.writelines(r1)
            o2.writelines(r2)
            n += 1
print(f"wrote {n} pairs to {out1}", file=sys.stderr)
