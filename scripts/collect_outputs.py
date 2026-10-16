#!/usr/bin/env python3
"""Copy selected files out of a finished Nextflow run's work directories into a flat folder.

usage: collect_outputs.py <mini|full> <tasks.tsv> <dest_dir>
  tasks.tsv = output of `nextflow log <run> -f name,workdir`

Every rule is (process-name regex, file-name regex, destination template).
The destination is formatted with {s} (task tag, i.e. the sample id or 'all'/'transcriptome'),
{f} (the source file name), {base} (same, kept for readability) and {stem} (file name without extension); `{rest}` is the source file name without a leading '<tag>_' or '<tag>.'.
Files that are staged inputs (symlinks created by Nextflow) are only taken if the rule says so via the
'in:' prefix in the file-name regex.
"""
import os
import re
import shutil
import sys

mode, tasks_file, dest = sys.argv[1:4]
os.makedirs(dest, exist_ok=True)

P = "NFCORE_CIRCRNA:CIRCRNA:"
RULES = {
    "mini": [
        # STAR 2-pass (BSJ_DETECTION:STAR2PASS only, not the circtools mate-wise STAR runs)
        (P + r"BSJ_DETECTION:STAR2PASS:PASS_2", r"\.Chimeric\.out\.junction$", "{s}.star.Chimeric.out.junction"),
        (P + r"BSJ_DETECTION:STAR2PASS:PASS_2", r"\.Chimeric\.out\.sam$", "{s}.star.Chimeric.out.sam"),
        (P + r"BSJ_DETECTION:STAR2PASS:PASS_2", r"\.SJ\.out\.tab$", "{s}.star.SJ.out.tab"),
        (P + r"BSJ_DETECTION:STAR2PASS:PASS_2", r"\.Aligned\.out\.bam$", "{s}.star.Aligned.out.bam"),
        # CIRI2
        (P + r"BSJ_DETECTION:CIRI:CIRI2", r"^{s}\.txt$", "{s}.ciri.txt"),
        (P + r"BSJ_DETECTION:CIRI:CIRI2", r"^{s}\.sam$", "{s}.ciri.sam"),
        (P + r"BSJ_DETECTION:CIRI:READLENGTH", r"_readlength\.txt$", "{s}.ciri.readlength.txt"),
        (P + r"BSJ_DETECTION:CIRI:FASTP", r"_R1\.fastp\.fastq\.gz$", "{s}.ciri.fixed_length_1.fastq.gz"),
        (P + r"BSJ_DETECTION:CIRI:FASTP", r"_R2\.fastp\.fastq\.gz$", "{s}.ciri.fixed_length_2.fastq.gz"),
        # psirc BSJ (output directory without the bulky mapping/ folder)
        (P + r"BSJ_DETECTION:PSIRC:BSJ", r"output/(candidate_circ_junctions\.(bed|fa)|candidate_circ_supporting_reads\.sam|junction_names_output_to_raw\.tsv)$", "{s}.psirc_bsj.{base}"),
        # BED files
        (P + r"BSJ_DETECTION:BED_ADD_SAMPLE_TOOL", r"\.filtered\.bed$", "{s}.{rest}"),
        (P + r"BSJ_DETECTION:BED_ADD_SAMPLE_TOOL", r"\.meta\.bed$", "{s}.{rest}"),
        (P + r"BSJ_DETECTION:COMBINEBEDS_READS", r"_reads\.bed$", "{s}.reads.bed"),
        (P + r"BSJ_DETECTION:COMBINE_TOOLS_PER_SAMPLE", r"\.combined\.bed$", "{f}"),
        (P + r"BSJ_DETECTION:COMBINE_SAMPLES", r"^all\.combined\.bed$", "{f}"),
        (P + r"BSJ_DETECTION:ANNOTATE_COMBINED:RENAME", r"\.annotated\.bed12$", "{f}"),
        (P + r"BSJ_DETECTION:ANNOTATE_COMBINED:BED2GTF", r"\.annotated\.gtf$", "{f}"),
        (P + r"BSJ_DETECTION:ANNOTATE_COMBINED:GET_FASTA", r"^all\.fa$", "{f}"),
        (P + r"BSJ_DETECTION:ANNOTATE_PER_SAMPLE:RENAME", r"\.annotated\.bed12$", "{f}"),
        (P + r"BSJ_DETECTION:ANNOTATE_PER_SAMPLE:BED2GTF", r"\.annotated\.gtf$", "{f}"),
        (P + r"BSJ_DETECTION:ANNOTATE_PER_SAMPLE:GET_FASTA", r"^{s}\.fa$", "{s}.annotated.fa"),
        # genome / transcriptome
        (P + r"PREPARE_GENOME:SAMTOOLS_FAIDX", r"\.fai$", "chrI_mini.fa.fai"),
        (P + r"COMBINE_TRANSCRIPTOMES:TRANSCRIPTOME", r"^transcriptome\.fasta$", "{f}"),
        (P + r"COMBINE_TRANSCRIPTOMES:EXCLUDE_OVERLONG_TRANSCRIPTS", r"^transcriptome\.filtered\.gtf$", "{f}"),
        # miRNA
        (P + r"MIRNA_PREDICTION:DESEQ2_NORMALIZATION", r"\.normalized_counts\.tsv$", "mirna_expression.normalized.tsv"),
        (P + r"MIRNA_PREDICTION:MIRNA_FILTERING", r"\.normalized_counts_filtered\.tsv$", "mirna_expression.filtered.tsv"),
        (P + r"MIRNA_PREDICTION:MIRNA_BINDINGSITES:MAJORITY_VOTE", r"\.majority\.tsv$", "mirna.binding_sites.majority.tsv"),
        (P + r"MIRNA_PREDICTION:MIRNA_BINDINGSITES:MAJORITY_VOTE", r"\.targets\.tsv$", "mirna.binding_sites.targets.tsv"),
        # quantification
        (P + r"QUANTIFICATION:PSIRC_QUANT:MERGE_EXPERIMENTS", r"^experiments\.merged\.rds$", "psirc.{f}"),
        (P + r"QUANTIFICATION:PSIRC_QUANT:JOIN_GENE_COUNTS", r"^gene_counts\.csv$", "psirc.{f}"),
        (P + r"QUANTIFICATION:PSIRC_QUANT:JOIN_GENE_TPM", r"^gene_tpm\.csv$", "psirc.{f}"),
        (P + r"QUANTIFICATION:PSIRC_QUANT:JOIN_TX_COUNTS", r"^tx_counts\.csv$", "psirc.{f}"),
        (P + r"QUANTIFICATION:PSIRC_QUANT:JOIN_TX_TPM", r"^tx_tpm\.csv$", "psirc.{f}"),
        (P + r"QUANTIFICATION:PSIRC_QUANT:SPLIT_TYPES_COUNTS", r"^(circular|linear)\.tsv$", "psirc.{stem}_tx_counts.tsv"),
        (P + r"QUANTIFICATION:CIRIQUANT:MAIN", r"/{s}\.gtf$", "{s}.ciriquant.gtf"),
        (P + r"QUANTIFICATION:CIRIQUANT:MAIN", r"/gene/{s}_out\.gtf$", "{s}.stringtie.gtf"),
        (P + r"QUANTIFICATION:CIRIQUANT:EXTRACT_CIRC", r"\.circ\.tsv$", "{s}.ciriquant.circ.tsv"),
        (P + r"QUANTIFICATION:CIRIQUANT:EXTRACT_GENES", r"\.gene\.tsv$", "{s}.ciriquant.gene.tsv"),
        (P + r"QUANTIFICATION:CIRIQUANT:JOIN_GENE", r"^gene\.joined\.tsv$", "ciriquant.gene_tpm.tsv"),
        (P + r"QUANTIFICATION:CIRIQUANT:JOIN_CIRC", r"^circ\.joined\.tsv$", "ciriquant.circ_cpm.tsv"),
    ],
    "full": [
        (P + r"QUANTIFICATION:PSIRC_QUANT:JOIN_GENE_COUNTS", r"^gene_counts\.csv$", "psirc.gene_counts.csv"),
        (P + r"QUANTIFICATION:PSIRC_QUANT:SPLIT_TYPES_COUNTS", r"^circular\.tsv$", "psirc.circular_tx_counts.tsv"),
        (P + r"QUANTIFICATION:CIRIQUANT:JOIN_GENE", r"^gene\.joined\.tsv$", "ciriquant.gene_tpm.tsv"),
        (P + r"QUANTIFICATION:CIRIQUANT:JOIN_CIRC", r"^circ\.joined\.tsv$", "ciriquant.circ_cpm.tsv"),
        (P + r"QUANTIFICATION:CIRIQUANT:MAIN", r"/{s}\.gtf$", "{s}.ciriquant.gtf"),
        (P + r"QUANTIFICATION:CIRIQUANT:MAIN", r"/gene/{s}_out\.gtf$", "{s}.stringtie.gtf"),
    ],
}

tasks = [l.rstrip("\n").split("\t") for l in open(tasks_file) if "\t" in l]
copied = []
for proc_re, file_re, dest_t in RULES[mode]:
    for name, workdir in tasks:
        if not re.search("^" + proc_re, name):
            continue
        m = re.search(r"\(([^)]*)\)\s*$", name)
        tag = m.group(1) if m else ""
        sample = tag
        # all files below the work dir (excluding Nextflow's own dot files and STAR scratch dirs)
        for root, dirs, files in os.walk(workdir):
            dirs[:] = [d for d in dirs if not d.startswith(".") and not d.endswith("_STARgenome") and d != "mapping"]
            for fn in files:
                if fn.startswith("."):
                    continue
                rel = os.path.relpath(os.path.join(root, fn), workdir)
                pat = file_re.replace("{s}", re.escape(sample))
                if not re.search(pat, rel if "/" in pat else fn):
                    continue
                rest = re.sub(r"^" + re.escape(sample) + r"[_.]", "", fn)
                base = fn
                out = dest_t.format(s=sample, f=fn, rest=rest, base=base, stem=os.path.splitext(fn)[0])
                src = os.path.join(root, fn)
                shutil.copyfile(src, os.path.join(dest, out))  # follows symlinks
                copied.append(out)
print("\n".join(sorted(set(copied))))
