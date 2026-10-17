# Test data for the nf-core/circrna module and subworkflow tests

Everything below lives in the fork `nictru/test-datasets` (branch `circrna`) and is read through
`params.pipelines_testdata_base_path` (set in `tests/nextflow.config` of nf-core/circrna to
`https://raw.githubusercontent.com/nictru/test-datasets/<commit-sha>/`). Always pin the commit SHA, never the branch.

In a test: `file(params.pipelines_testdata_base_path + <path from the tables below>, checkIfExists: true)`, e.g. the path `mini/reference/chrI_mini.fa`.

The upstream data (`fastq/`, `reference/`, `bwa/`, `hisat2/`, `miranda/`, `circexplorer2/`, `samples.csv`, `phenotype.csv`)
stays untouched at the root of the branch; the user-facing `-profile test` keeps reading **upstream** (see "Base path" below).

## Overview

| Set | What | Size |
| --- | --- | --- |
| mini | 800 kb region of C. elegans chrI, 2 samples x 20 000 read pairs, frozen outputs of one run with all tools | ~48 MB |
| full stats | quantification outputs of the 6-sample `-profile test` data | ~29 MB |

**Region**: `chrI:1125001-1925000` (800 kb). It is the window with the most circexplorer2 BSJs (427 unique BSJs from the 6
samples) that lie completely inside genes that are completely inside the window. The window is cut out of `reference/chrI.fa`
and `reference/chrI.gtf` and **rebased to start at 1**; the sequence is still called `chrI` (so a BSJ at `chrI:a-b` in the
mini data is at `chrI:a+1125000-b+1125000` in upstream coordinates). 121 genes / 2293 GTF features.

**Samples** (paired-end, unstranded): `cel_N2_1` (from upstream `N2_rep1`) and `fust1_1` (from `fust1_rep1`).
Meta used by the pipeline: `[id: 'cel_N2_1', single_end: false]`, strandedness `auto`.
Reads: all read pairs that STAR reports as chimeric (BSJ) reads with both breakpoints in the region (387 / 381 pairs) plus a
deterministic md5-ordered subsample of the pairs whose mates both align in the region, 20 000 pairs per sample in total.

## Base path

`conf/test.config` builds its URLs from `params.pipelines_testdata_base_path` while the profile is evaluated, i.e. with the
default of `nextflow.config` (`https://raw.githubusercontent.com/nf-core/test-datasets/`), **not** with the value from
`tests/nextflow.config` (which is applied afterwards via `-c`). Checked with `nextflow -c <config with the fork path> config -profile test .`:
`input`, `fasta`, `gtf`, `mature`, `phenotype` still resolve to `.../nf-core/test-datasets/circrna/...`. So the fork path only
affects module/subworkflow tests (which read `params.pipelines_testdata_base_path` inside the test), the pipeline test keeps using upstream data.

## Reference, reads and miRNA inputs

@@REFERENCE_TABLE@@

`mini/mirna_expression.tsv` has the columns `miRNA, cel_N2_1, cel_N2_2, cel_N2_3, fust1_1, fust1_2, fust1_3` (random counts, seed 42).
`MIRNA_FILTERING` refuses fewer than 5 sample columns (and reads the counts from column 5 on), so the table cannot be restricted to the 2 mini samples.

`mini/samples_mini.csv` points at raw.githubusercontent.com URLs of the fork that use the branch name `circrna` (because the SHA is only known after the commit).
Tests should build their input channels from `params.pipelines_testdata_base_path` and the fastq paths in the first table instead.

## Frozen intermediates

One run of the pipeline on the mini set (generator step "mini run A/B") with `--tools ciri,circexplorer2,find_circ,circrna_finder,mapsplice,circtools,segemehl,psirc`,
`--fli_tools psirc,circtools,cirifull`, `--quantification_tools psirc,ciriquant,sum,max --bootstrap_samples 3`, `--mirna_tools miranda,targetscan`, `--min_tools 1`.
File names are `<sample>.<what>.<ext>` (or `all.*`, `psirc.*`, `ciriquant.*`, `mirna*.*` for files that are not per sample).
Channel names in the "What" column are the emits of `BSJ_DETECTION`, `COMBINE_TRANSCRIPTOMES` and `QUANTIFICATION`.

@@INTERMEDIATES_TABLE@@

Meta fields the subworkflows expect around these files:

- per sample files: `[id: '<sample>', single_end: false]`; `bed_per_sample_tool`: `[id: '<sample>', tool: '<tool>']` (+ the BED has the sample/tool columns);
  the fixed-length reads: `[id: '<sample>', target_length: <int from *.ciri.readlength.txt>]`;
- combined files: `[id: 'all']`; transcriptome files: `[id: 'transcriptome']`;
- psirc quantification matrices: `[id: 'gene_counts' | 'tx_counts' ..., quantification: 'psirc']`; CIRIquant matrices `[id: 'gene' | 'circ']` (see "Findings" for why the metas are inconsistent).

## Full stats (inputs for `STATISTICAL_TESTS`)

Produced by the 6-sample `-profile test` data (`--tools circexplorer2 --quantification_tools ciriquant,psirc --bootstrap_samples 3`,
`mature = null`, upstream `bwa/` and `hisat2/` indices passed as `--bwa/--hisat2`, see "Findings" 2).
Channel shapes of `STATISTICAL_TESTS`:

- `ch_gene_counts` / `ch_circ_counts`: `[meta, file]` with `meta.quantification` in `psirc`, `ciriquant` (psirc: `psirc.gene_counts.csv` / `psirc.circular_tx_counts.tsv`; ciriquant: `ciriquant.gene_tpm.tsv` / `ciriquant.circ_cpm.tsv`). **Both channels are joined on the whole meta, so a test has to give gene and circ the same meta** (e.g. `[id: 'counts', quantification: 'psirc']`);
- `ch_ciriquant` / `ch_stringtie`: `[[id: '<sample>'], gtf]` for the 6 samples `cel_N2_1..3`, `fust1_1..3`;
- `ch_phenotype`: `[[id: 'phenotype'], phenotype.csv]`.

@@STATS_TABLE@@

## Tools that are missing or degraded

State of the tools when the data were generated. The `FLI_DETECTION`, `CIRIVIS_UNIFY` and `jccirc` problems listed here are fixed in nf-core/circrna PR #242, but no frozen FLI outputs are shipped: the cirifull module tests start from the CIRI inputs listed above.

| Tool / step | State on the mini set | Consequence |
| --- | --- | --- |
| BSJ: circexplorer2, circrna_finder, circtools, ciri, find_circ, mapsplice, segemehl | worked, intermediates present | real tests possible |
| BSJ: psirc | ran, found **0** junctions (also with `-m 5 -u 100`, other library types) | `*.psirc_bsj.*` files are empty; use `-stub` for psirc BSJ/FLI-from-BSJ |
| FLI: psirc, circtools | worked inside the run, but `FLI_DETECTION` mixes the psirc FASTA into the genome FASTA that circtools gets (Finding 1) | no frozen FLI outputs shipped |
| FLI: cirifull | `CIRIAS`/`RO2`/`MERGE`/`BUILD_LIST`/`CIRI_VIS` never start because `FLI_DETECTION` overwrites its `ch_fasta` input with an empty channel (Finding 1, which hides a second bug in the `CIRIVIS_UNIFY` output); the CIRI inputs (`*.ciri.txt`, `*.ciri.sam`, fixed-length reads) are present | the cirifull module tests start directly from these files |
| FLI: jccirc | not run: `jccirc` is not allowed by the `fli_tools` pattern of `nextflow_schema.json` | stub only |
| Quantification: CIRIquant | runs, but quantifies **no** circRNA on the mini set (header-only GTF/matrix; 5 circRNAs on the 6-sample data) | module test passes on empty output; the 6-sample files in the stats set are non-trivial |
| miRNA: miranda, targetscan | worked | |
| miRNA: expression branch (`mirna_expression`) | `COMPUTE_CORRELATIONS` fails with "No transcripts left after filtering" (fishpond needs >= 3 samples) | `DESEQ2_NORMALIZATION` and `MIRNA_FILTERING` outputs are shipped; correlations: `-stub` + TODO |

## Findings while generating the data

Behaviour of the pipeline at the time the data were generated. Findings 1, 4 and 5 (and the `jccirc` schema entry above) are fixed in nf-core/circrna PR #242; 2, 3 and 6 are tracked as issues in nf-core/circrna.

1. `subworkflows/local/fli_detection/main.nf`: the local `ch_fasta = channel.empty()` shadows the input `ch_fasta` (genome FASTA) and is later mixed with the FLI FASTAs, so CIRIFULL/JCCIRC get an empty channel (cirifull silently stops after `RO1`) and CIRCTOOLS gets the wrong FASTA. Renaming the accumulator (e.g. `ch_fli_fasta`) fixes it; with that patch cirifull runs through `RO2`, `MERGE`, `BUILD_LIST` and `CIRI_VIS` on the mini data. It then fails in the `CIRIFULL` subworkflow because `CIRIVIS_UNIFY` emits `bed` as a bare `path` (no `[meta, path]` tuple) while the subworkflow does `.map { meta, bed -> ... }` ("Invalid method invocation `call`"): the module needs `tuple val(meta), path("${prefix}.bed")`.
2. `PREPARE_GENOME` only builds the bwa/hisat2 indices when `ciri` is a detection tool (or `--bwa/--hisat2` are given). With `--tools circexplorer2 --quantification_tools ciriquant` (the `-profile test` default) `CIRIQUANT` silently never runs because its index channels are empty.
3. `STATISTICAL_TESTS` joins `ch_gene_counts` and `ch_circ_counts` on the full meta, but the quantification subworkflows emit different metas (`gene_counts` vs `tx_counts`, `gene` vs `circ`), so `CIRCTEST_PREPARE`/`CIRCTEST_CIRCTEST` never run in the pipeline (only `CIRIQUANT_DE` does).
4. `BSJ_DETECTION` fails at the `emit:` block (`STAR2PASS.out` undefined) when no STAR-based tool is selected (e.g. `--tools ciri`).
5. `MIRNA_BINDINGSITES` fails with `--mirna_tools miranda` alone ("Access to 'TARGETSCAN.out' is undefined", `mirna_bindingsites/main.nf` line 66); the tests must select both tools or the subworkflow needs a guard.
6. `MIRNA_FILTERING` needs >= 5 sample columns and loops from column 5; `COMPUTE_CORRELATIONS` still reads `params.mirna_correlation` inside its template and needs >= 3 samples.

## Reproducing everything

```bash
git clone -b circrna git@github.com:nictru/test-datasets.git && cd test-datasets
# needs: nextflow >= 25.10, apptainer (or PROFILE=docker), samtools, python3; ~10 GB scratch, 16+ cpus
PROFILE=apptainer bash scripts/generate_testdata.sh /path/to/nf-core/circrna /path/to/scratch
python3 scripts/write_manifest.py . README.md <commit-sha>
```

`scripts/generate_testdata.sh` runs: (1) the 6-sample test profile -> `full stats`; (2) `select_region.py` (window with most BSJs);
(3) `make_mini_reference.py`, `pick_read_pairs.py` (mini reference and fastqs, byte-reproducible given the same STAR BAMs; STAR multithreading can change a handful of pairs between runs);
(4) mini run A with all tools and the complete `reference/mature.fa`; (5) `select_mirnas.py` (20 miRNAs, expression table); (6) mini run B (20 miRNAs) and run C (+ expression table, expected to fail in `COMPUTE_CORRELATIONS`);
(7) `collect_outputs.py` copies the intermediates out of the Nextflow work dirs (`nextflow log -f name,workdir`).
Pipeline revision used: nf-core/circrna `dev` after "refactor: move local subworkflows into folders and make them params-free" (params-free subworkflows).
Run times (16 cpus, 100 GB, apptainer, warm image cache): 6-sample run ~7 min, mini run A ~6 min, mini runs B/C ~3 min each.
