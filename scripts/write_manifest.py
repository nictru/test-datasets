#!/usr/bin/env python3
"""Write the test data manifest into README.md (between the TESTDATA markers) from the files that exist in the checkout.

usage: write_manifest.py <testdata_dir> <README.md> [<sha-or-branch>]
"""
import os
import re
import sys

td, out = sys.argv[1:3]
ref = sys.argv[3] if len(sys.argv) > 3 else "circrna"

def files(sub):
    root = os.path.join(td, sub)
    return sorted(os.path.relpath(os.path.join(r, f), td) for r, _, fs in os.walk(root) for f in fs)

def size(p):
    n = os.path.getsize(os.path.join(td, p))
    return f"{n / 1e6:.1f} MB" if n >= 1e6 else f"{max(1, n // 1000)} kB" if n >= 1000 else f"{n} B"

# (regex on the file name, what it is, producing process, consuming test packages)
INTER = [
    (r"\.star\.Chimeric\.out\.junction$", "STAR 2-pass `Chimeric.out.junction` (BSJ evidence)", "`STAR2PASS:PASS_2`"),
    (r"\.star\.Chimeric\.out\.sam$", "STAR 2-pass `Chimeric.out.sam` (written because `circrna_finder` is selected)", "`STAR2PASS:PASS_2`"),
    (r"\.star\.SJ\.out\.tab$", "STAR 2-pass splice junctions", "`STAR2PASS:PASS_2`"),
    (r"\.star\.Aligned\.out\.bam$", "STAR 2-pass unsorted BAM (`ch_star_bam`)", "`STAR2PASS:PASS_2`"),
    (r"\.ciri2\.txt$", "CIRI2 result table (`ciri2_txt`)", "`CIRI2:DETECT`"),
    (r"\.ciri2\.sam$", "BWA-MEM alignment fed to CIRI2 (`ciri2_sam`)", "`CIRI_BWA_MEM`"),
    (r"\.ciri2\.readlength\.txt$", "5th-percentile read length; becomes `meta.target_length` of the fixed-length reads", "`CIRI_READS:READLENGTH`"),
    (r"\.ciri2\.fixed_length_[12]\.fastq\.gz$", "reads trimmed to equal length (`reads_fixed_length`, meta `[id, single_end:false, target_length:<readlength>]`)", "`CIRI_READS:FASTP`"),
    (r"\.psirc_bsj\.", "psirc BSJ output directory content (`output/`, without `mapping/`); **empty: psirc finds 0 junctions**", "`PSIRC:BSJ`"),
    (r"\.[a-z0-9_]+\.filtered\.bed$", "per-sample, per-tool unified BED after the `bsj_reads` filter (input of `BED_ADD_SAMPLE_TOOL`)", "`FILTER_BSJS`"),
    (r"\.[a-z0-9_]+\.meta\.bed$", "per-sample, per-tool BED with sample/tool columns (`bed_per_sample_tool`, meta `[id:<sample>, tool:<tool>]`)", "`BED_ADD_SAMPLE_TOOL`"),
    (r"\.reads\.bed$", "read-level agreement BED (`bed_reads`, meta `[id:<sample>]`); nearly empty on the mini data", "`COMBINEBEDS_READS`"),
    (r"^(cel_N2_1|fust1_1)\.combined\.bed$", "per-sample combined BED (`bed_per_sample`)", "`COMBINE_TOOLS_PER_SAMPLE`"),
    (r"^all\.combined\.bed$", "all samples/tools combined BED (`bed`)", "`COMBINE_SAMPLES`"),
    (r"\.annotated\.bed12$", "annotated BED12 (`bed12`)", "`ANNOTATE_*:RENAME`"),
    (r"\.annotated\.gtf$", "annotated circRNA GTF (`gtf`)", "`ANNOTATE_*:BED2GTF`"),
    (r"^all\.fa$", "FASTA of the combined circRNAs (`fasta`)", "`ANNOTATE_COMBINED:GET_FASTA`"),
    (r"\.annotated\.fa$", "FASTA of the per-sample circRNAs", "`ANNOTATE_PER_SAMPLE:GET_FASTA`"),
    (r"^chrI_mini\.fa\.fai$", "`samtools faidx` index of the mini genome (`faidx`)", "`PREPARE_GENOME:SAMTOOLS_FAIDX`"),
    (r"^transcriptome\.fasta$", "combined linear + circular transcriptome FASTA", "`COMBINE_TRANSCRIPTOMES:TRANSCRIPTOME`"),
    (r"^transcriptome\.filtered\.gtf$", "combined transcriptome GTF (overlong transcripts removed)", "`COMBINE_TRANSCRIPTOMES:EXCLUDE_OVERLONG_TRANSCRIPTS`"),
    (r"^psirc\.experiments\.merged\.rds$", "psirc-quant merged `SummarizedExperiment` (`rds`)", "`PSIRC_QUANT:MERGE_EXPERIMENTS`"),
    (r"^psirc\..*\.(csv|tsv)$", "psirc-quant count/TPM matrices (gene, transcript, circular, linear)", "`PSIRC_QUANT:JOIN_*` / `SPLIT_TYPES_COUNTS`"),
    (r"\.ciriquant\.gtf$", "CIRIquant per-sample circRNA GTF (`raw`); header only on the mini data", "`CIRIQUANT:MAIN`"),
    (r"\.stringtie\.gtf$", "CIRIquant per-sample StringTie GTF (`stringtie`)", "`CIRIQUANT:MAIN`"),
    (r"\.ciriquant\.circ\.tsv$", "CIRIquant per-sample circ table", "`CIRIQUANT:EXTRACT_CIRC`"),
    (r"\.ciriquant\.gene\.tsv$", "CIRIquant per-sample gene table", "`CIRIQUANT:EXTRACT_GENES`"),
    (r"^ciriquant\.(gene_tpm|circ_cpm)\.tsv$", "CIRIquant joined matrix (circ matrix is header only on the mini data)", "`CIRIQUANT:JOIN_GENE` / `JOIN_CIRC`"),
    (r"^mirna_expression\.normalized\.tsv$", "DESeq2-normalised miRNA expression (6 sample columns)", "`DESEQ2_NORMALIZATION`"),
    (r"^mirna_expression\.filtered\.tsv$", "filtered normalised miRNA expression", "`MIRNA_FILTERING`"),
    (r"^mirna\.binding_sites\.majority\.tsv$", "miRNA -> circRNA binding sites after majority vote (`min_tools` 1)", "`MAJORITY_VOTE`"),
    (r"^mirna\.binding_sites\.targets\.tsv$", "miRNA -> target list (`MIRNA_TARGETS` format)", "`MAJORITY_VOTE`"),
]

def describe(name):
    base = os.path.basename(name)
    for rx, what, proc in INTER:
        if re.search(rx, base):
            return what, proc
    return "(undocumented)", ""

def table(rows):
    lines = ["| File | Size | What | Produced by |", "| --- | --- | --- | --- |"]
    for p, what, proc, *_ in rows:
        lines.append(f"| `{p}` | {size(p)} | {what} | {proc} |")
    return "\n".join(lines)

ref_files = [
    ("mini/reference/chrI_mini.fa", "mini genome: one `chrI` sequence = chrI region rebased to 1", "`scripts/make_mini_reference.py`"),
    ("mini/reference/chrI_mini.gtf", "mini annotation: genes fully inside the region, coordinates rebased", "`scripts/make_mini_reference.py`"),
    ("mini/reference/mature_mini.fa", "20 miRNAs with the most predicted sites on the mini circRNAs (miranda AND targetscan)", "`scripts/select_mirnas.py`"),
    ("mini/mirna_expression.tsv", "synthetic miRNA counts, 6 sample columns (see below)", "`scripts/select_mirnas.py`"),
    ("mini/phenotype.csv", "`sample,condition` for the two mini samples", "`scripts/generate_testdata.sh`"),
    ("mini/samples_mini.csv", "samplesheet whose URLs point at the `circrna` branch of the fork", "`scripts/generate_testdata.sh`"),
    ("mini/fastq/N2_rep1_1.fastq.gz", "sample `cel_N2_1`, mate 1", "`scripts/pick_read_pairs.py`"),
    ("mini/fastq/N2_rep1_2.fastq.gz", "sample `cel_N2_1`, mate 2", "`scripts/pick_read_pairs.py`"),
    ("mini/fastq/fust1_rep1_1.fastq.gz", "sample `fust1_1`, mate 1", "`scripts/pick_read_pairs.py`"),
    ("mini/fastq/fust1_rep1_2.fastq.gz", "sample `fust1_1`, mate 2", "`scripts/pick_read_pairs.py`"),
]
inter_rows = [(p, *describe(p)) for p in files("mini/intermediates")]
stats_desc = [
    (r"psirc\.gene_counts\.csv", "psirc gene counts, `meta.quantification = 'psirc'` (columns: tx/gene id, 6 samples)"),
    (r"psirc\.circular_tx_counts\.tsv", "psirc circular transcript counts (`circ`), `meta.quantification = 'psirc'`"),
    (r"ciriquant\.gene_tpm\.tsv", "CIRIquant gene TPM matrix (`gene`)"),
    (r"ciriquant\.circ_cpm\.tsv", "CIRIquant circRNA CPM matrix (`circ`; 5 circRNAs)"),
    (r"\.ciriquant\.gtf$", "CIRIquant per-sample circRNA GTF (`ciriquant`, meta `[id:<sample>]`)"),
    (r"\.stringtie\.gtf$", "CIRIquant per-sample StringTie GTF (`stringtie`, meta `[id:<sample>]`)"),
    (r"phenotype\.csv", "`sample,condition` of the 6 samples (3x N2, 3x fust1)"),
]
def sdesc(p):
    for rx, d in stats_desc:
        if re.search(rx, p):
            return d
    return ""
stats_rows = [(p, sdesc(p), "full test-profile run" if "phenotype" not in p else "copy of the upstream file") for p in files("full/stats")]

tmpl = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "manifest_template.md")).read()
doc = (tmpl.replace("@@REF@@", ref)
       .replace("@@REFERENCE_TABLE@@", table(ref_files))
       .replace("@@INTERMEDIATES_TABLE@@", table(inter_rows))
       .replace("@@STATS_TABLE@@", table(stats_rows)))
BEGIN, END = "<!-- BEGIN TESTDATA (generated by scripts/write_manifest.py) -->", "<!-- END TESTDATA -->"
lines, fenced = [], False
for line in doc.split("\n"):
    if line.startswith("```"):
        fenced = not fenced
    lines.append("##" + line if line.startswith("#") and not fenced else line)
body = "\n".join(lines)
readme = open(out).read()
block = f"{BEGIN}\n{body.rstrip()}\n{END}"
if BEGIN in readme:
    readme = re.sub(re.escape(BEGIN) + r".*?" + re.escape(END), lambda _: block, readme, flags=re.S)
else:
    readme = readme.rstrip() + "\n\n" + block + "\n"
open(out, "w").write(readme)
print(f"{out}: {len(ref_files) + len(inter_rows) + len(stats_rows)} files documented", file=sys.stderr)
undoc = [p for p, w, *_ in inter_rows if w == "(undocumented)"]
if undoc:
    print("UNDOCUMENTED:", undoc, file=sys.stderr)
