#!/usr/bin/env bash
# Generate the nf-core/circrna module/subworkflow test data (mini/** and full/stats/**).
#
# usage: bash scripts/generate_testdata.sh <pipeline_dir> <work_dir>
#
#   <pipeline_dir>  checkout of nf-core/circrna (dev, with params-free subworkflows)
#   <work_dir>      scratch directory for Nextflow work dirs and run outputs (~10 GB)
#
# Writes into the checkout this script lives in (override with TESTDATA_DIR):
#   mini/reference/{chrI_mini.fa,chrI_mini.gtf,mature_mini.fa}, mini/fastq/*, mini/mirna_expression.tsv,
#   mini/samples_mini.csv, mini/phenotype.csv, mini/intermediates/*, full/stats/*
#
# Requirements: nextflow, an apptainer/singularity/docker runtime (PROFILE), samtools, python3.
# Environment: PROFILE (default apptainer), REGION_SIZE (default 800000), MAX_PAIRS (default 20000),
#              CPUS (default 16), MEM (default 100.GB), BRANCH_URL (raw URL prefix used in samples_mini.csv)
set -euo pipefail

PIPE=$(realpath "${1:?usage: generate_testdata.sh <pipeline_dir> <work_dir>}")
WORK=$(realpath -m "${2:?usage: generate_testdata.sh <pipeline_dir> <work_dir>}")
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TD=${TESTDATA_DIR:-$(dirname "$HERE")}
PROFILE=${PROFILE:-apptainer}
REGION_SIZE=${REGION_SIZE:-800000}
MAX_PAIRS=${MAX_PAIRS:-20000}
CPUS=${CPUS:-16}
MEM=${MEM:-100.GB}
BRANCH_URL=${BRANCH_URL:-https://raw.githubusercontent.com/nictru/test-datasets/circrna/}
CHROM=chrI

mkdir -p "$WORK/full" "$WORK/mini" "$TD/mini/reference" "$TD/mini/fastq" "$TD/mini/intermediates" "$TD/full/stats"

cat > "$WORK/resources.config" <<EOF
process.resourceLimits = [cpus: $CPUS, memory: '$MEM', time: '12.h']
EOF
cat > "$WORK/full/no_mature.config" <<EOF
params.mature = null
EOF

# Run nextflow; image pulls occasionally fail when several pulls race, so retry with -resume.
nf() {
    local prefix=$1 n=0
    shift
    until nextflow run "$@" -name "${prefix}_$(date +%s)"; do
        n=$((n + 1)); [ $n -ge 3 ] && return 1
    done
}

# ---------------------------------------------------------------------------------------------
# 1. Full 6-sample chrI run (the `test` profile): source of the region choice and of full/stats
#    quantification_tools ciriquant only runs when the bwa/hisat2 indices are provided or `ciri`
#    is a detection tool, so the upstream prebuilt indices are passed in.
# ---------------------------------------------------------------------------------------------
cd "$WORK/full"
nf full "$PIPE" -profile "test,$PROFILE" -resume \
    -c "$WORK/resources.config" -c no_mature.config -w "$WORK/full/work" \
    --outdir "$WORK/full/out" --tools circexplorer2 --quantification_tools ciriquant,psirc \
    --bootstrap_samples 3 --bwa "$TD/bwa" --hisat2 "$TD/hisat2" --save_intermediates
FULL_RUN=$(nextflow log -q | tail -1)
nextflow log "$FULL_RUN" -f name,workdir > tasks.tsv
rm -rf "$WORK/full/stats" && python3 "$HERE/collect_outputs.py" full tasks.tsv "$WORK/full/stats" > /dev/null
cp "$WORK"/full/stats/* "$TD/full/stats/"
cp "$TD/phenotype.csv" "$TD/full/stats/phenotype.csv"

# ---------------------------------------------------------------------------------------------
# 2. Region choice: window with most BSJs (circexplorer2, all 6 samples) inside complete genes
# ---------------------------------------------------------------------------------------------
cat "$WORK"/full/out/3_bsj_detection/tools/circexplorer2/unified/*.bed > "$WORK/full/all_circexplorer2.bed"
read -r _ START END NBSJ < <(python3 "$HERE/select_region.py" "$WORK/full/all_circexplorer2.bed" "$TD/reference/chrI.gtf" $CHROM "$REGION_SIZE" 5000)
echo "region: $CHROM:$START-$END ($NBSJ unique BSJs in complete genes)"

python3 "$HERE/make_mini_reference.py" "$TD/reference/chrI.fa" "$TD/reference/chrI.gtf" $CHROM "$START" "$END" \
    "$TD/mini/reference/chrI_mini.fa" "$TD/mini/reference/chrI_mini.gtf"

STAR_DIR="$WORK/full/out/3_bsj_detection/tools/star/2nd_pass"
for pair in cel_N2_1:N2_rep1 fust1_1:fust1_rep1; do
    s=${pair%%:*}; f=${pair##*:}
    python3 "$HERE/pick_read_pairs.py" "$STAR_DIR/$s/$s.Aligned.out.bam" "$STAR_DIR/$s/$s.Chimeric.out.junction" \
        $CHROM "$START" "$END" "$TD/fastq/${f}_1.fastq.gz" "$TD/fastq/${f}_2.fastq.gz" \
        "$TD/mini/fastq/${f}_1.fastq.gz" "$TD/mini/fastq/${f}_2.fastq.gz" "$MAX_PAIRS"
done

printf 'sample,condition\ncel_N2_1,N2\nfust1_1,fust1\n' > "$TD/mini/phenotype.csv"
{
    echo "sample,fastq_1,fastq_2"
    echo "cel_N2_1,${BRANCH_URL}mini/fastq/N2_rep1_1.fastq.gz,${BRANCH_URL}mini/fastq/N2_rep1_2.fastq.gz"
    echo "fust1_1,${BRANCH_URL}mini/fastq/fust1_rep1_1.fastq.gz,${BRANCH_URL}mini/fastq/fust1_rep1_2.fastq.gz"
} > "$TD/mini/samples_mini.csv"
{
    echo "sample,fastq_1,fastq_2"
    echo "cel_N2_1,$TD/mini/fastq/N2_rep1_1.fastq.gz,$TD/mini/fastq/N2_rep1_2.fastq.gz"
    echo "fust1_1,$TD/mini/fastq/fust1_rep1_1.fastq.gz,$TD/mini/fastq/fust1_rep1_2.fastq.gz"
} > "$WORK/mini/samples_local.csv"

# ---------------------------------------------------------------------------------------------
# 3. Mini run A: all BSJ/FLI/quantification tools + miRNA tools with the complete mature.fa
#    (only used to find miRNAs with binding sites on the mini circRNAs).
#    jccirc is not accepted by the schema pattern of --fli_tools, so it is not part of the run.
# ---------------------------------------------------------------------------------------------
cd "$WORK/mini"
MINI_ARGS=(--input "$WORK/mini/samples_local.csv" --fasta "$TD/mini/reference/chrI_mini.fa" --gtf "$TD/mini/reference/chrI_mini.gtf"
    --tools ciri,circexplorer2,find_circ,circrna_finder,mapsplice,circtools,segemehl,psirc
    --fli_tools psirc,circtools,cirifull --quantification_tools psirc,ciriquant,sum,max --bootstrap_samples 3
    --mirna_tools miranda,targetscan --min_tools 1 --save_intermediates --save_reference)
nf miniA "$PIPE" -profile "$PROFILE" -resume -c "$WORK/resources.config" -w "$WORK/mini/work" \
    --outdir "$WORK/mini/outA" "${MINI_ARGS[@]}" --mature "$TD/reference/mature.fa"

python3 "$HERE/select_mirnas.py" "$WORK/mini/outA/6_mirna_prediction/binding_sites" "$TD/reference/mature.fa" 20 \
    "$TD/mini/reference/mature_mini.fa" "$TD/mini/mirna_expression.tsv" \
    cel_N2_1 cel_N2_2 cel_N2_3 fust1_1 fust1_2 fust1_3

# ---------------------------------------------------------------------------------------------
# 4. Mini run B: with the 20 miRNAs (no expression table; completes)
#    Mini run C: additionally with the synthetic expression table. This run is EXPECTED TO FAIL in
#    COMPUTE_CORRELATIONS ("No transcripts left after filtering": fishpond needs >= 3 samples, the mini set
#    has 2). It is only used to obtain the DESEQ2_NORMALIZATION / MIRNA_FILTERING outputs.
#    MIRNA_FILTERING needs >= 5 sample columns, hence the table has the 6 samples of the full data set.
# ---------------------------------------------------------------------------------------------
nf miniB "$PIPE" -profile "$PROFILE" -resume -c "$WORK/resources.config" -w "$WORK/mini/work" \
    --outdir "$WORK/mini/outB" "${MINI_ARGS[@]}" --mature "$TD/mini/reference/mature_mini.fa"
nextflow log -q | tail -1 | xargs -I{} nextflow log {} -f name,workdir > tasks.tsv
nextflow run "$PIPE" -profile "$PROFILE" -resume -c "$WORK/resources.config" -w "$WORK/mini/work" -name "miniC_$(date +%s)" \
    --outdir "$WORK/mini/outC" "${MINI_ARGS[@]}" --mature "$TD/mini/reference/mature_mini.fa" \
    --mirna_expression "$TD/mini/mirna_expression.tsv" || echo "miniC failed as expected (COMPUTE_CORRELATIONS)"
nextflow log -q | tail -1 | xargs -I{} nextflow log {} -f name,workdir >> tasks.tsv
sort -u tasks.tsv -o tasks.tsv
rm -rf "$WORK/mini/intermediates" && python3 "$HERE/collect_outputs.py" mini tasks.tsv "$WORK/mini/intermediates" > "$WORK/mini/intermediates.list"
cp "$WORK"/mini/intermediates/* "$TD/mini/intermediates/"
echo "done: $(wc -l < "$WORK/mini/intermediates.list") intermediate files"
