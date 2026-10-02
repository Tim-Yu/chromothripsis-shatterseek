# chromothripsis-shatterseek

Chromothripsis calling with [ShatterSeek](https://github.com/parklab/ShatterSeek) (Cortés-Ciriano et al. 2020,
Nat Genet) for the GEL sarcoma SV / copy-number data, plus an improved version of the group's legacy
PCF-density chromothripsis caller.

Two input layouts are supported by two runner scripts that share one library and write identical output tables:

| Layout | Script | Input |
|---|---|---|
| Per-sample files (circos pipeline, `build_sv_cnv_matrix.R` / `circos_plot_new.R`) | `scripts/run_shatterseek_persample.R` | consensus bedpe `<SID>.somatic.sv.bedpe` + Battenberg `<SID>_copynumber_1_gapfilled_LogR_new.txt` (or DRAGEN-style CN export) |
| Cohort tables (legacy `chromothripsis.r`) | `scripts/run_shatterseek_cohort.R` | `SVs.txt` (samplename chr1 strand1 pos1 chr2 strand2 pos2 svclass) + `CNAs.txt` (samplename chr start end nMaj1 nMin1 frac1 nMaj2 nMin2 frac2 SD ploidy) |

```
scripts/
  shatterseek_lib.R              shared: readers, harmonisation, ShatterSeek run, calling rules, plots
  run_shatterseek_persample.R    per-sample bedpe + CN files (single sample or --matrix batch)
  run_shatterseek_cohort.R       SVs.txt / CNAs.txt cohort tables
  run_gel_example.sh             GEL RE: patient list -> matrix -> ShatterSeek (per-sample runner)
  run_gel_legacy.sh              GEL RE: patient list -> matrix -> cohort tables -> improved legacy caller
  build_sv_cnv_matrix.R          patient list + SV/CN roots -> matrix (patientID sampleID SV_file_dir CN_file_dir)
  build_cohort_tables.R          matrix -> SVs.txt + CNAs.txt (legacy cohort format) from bedpe + Battenberg files
  gh_raw_fetch.sh                download files or whole public GitHub repos via github.com /raw/ addresses (no git, no API)
  install_shatterseek_RE.R       install ShatterSeek from the fetched source inside the research environment
legacy/
  chromothripsis_improved.r          improved legacy caller (same inputs/columns as the original + new statistics, tiered calls)
  chromothripsis_original_runnable.r original algorithm, unchanged, as a CLI (for comparison)
  reclassify_pui_table.r             apply the improved final rule to an existing legacy output table
  README_legacy.md                   algorithm of the original, why it missed events, what was changed
```

## 1. Install

R >= 4.0. ShatterSeek is not on CRAN.

```bash
# online
Rscript -e 'remotes::install_github("parklab/ShatterSeek", upgrade="never")'

# GEL research environment (no CRAN/GitHub; dependencies live in the shared library trees)
bash scripts/gh_raw_fetch.sh parklab/ShatterSeek -o ShatterSeek        # source via github.com /raw/ only
Rscript scripts/install_shatterseek_RE.R --src ShatterSeek             # installs into your personal library
#   the installer appends /tools/aws-workspace-ubuntu-apps/ce/R/4.5.3 and /tools/aws-workspace-apps/ce/R/4.2.1/
#   to .libPaths() (equivalent to the usual .libPaths(c(.libPaths(), "/tools/...")) lines) so the dependencies
#   are found; change them with --extra-libs dir1,dir2 and the target library with --lib <dir>
```
Dependencies: `BiocGenerics`, `graph`, `S4Vectors`, `GenomicRanges`, `IRanges`, `MASS`, `ggplot2`, `gridExtra`,
`foreach` (ShatterSeek); `copynumber`, `MASS`, `GenomicRanges` (legacy scripts). All are in the GEL shared trees.
The runner and legacy scripts set the library path automatically when the directories exist, equivalent to
```r
.libPaths(c("/home/byu/R/x86_64-pc-linux-gnu-library/4.5", .libPaths(),
            "/tools/aws-workspace-ubuntu-apps/ce/R/4.5.3", "/tools/aws-workspace-apps/ce/R/4.2.1/"))
```
(personal library with ShatterSeek first, shared trees with the dependencies last). Another personal library can be
given with `SHATTERSEEK_LIB=<dir>`, extra shared trees with `SHATTERSEEK_EXTRA_LIBS="dir1:dir2"`. Use the same
line in an interactive R session before `library(ShatterSeek)`.

```bash
git clone https://github.com/Tim-Yu/chromothripsis-shatterseek.git
```

**Restricted environments** (no git, no api.github.com; only github.com reachable): save `scripts/gh_raw_fetch.sh`
once (`curl -L -o gh_raw_fetch.sh https://github.com/Tim-Yu/chromothripsis-shatterseek/raw/main/scripts/gh_raw_fetch.sh`), then

```bash
# whole repository: the file list is read from the repo's tree pages, every file from its /raw/ address
bash gh_raw_fetch.sh Tim-Yu/chromothripsis-shatterseek -o chromothripsis-shatterseek
bash gh_raw_fetch.sh parklab/ShatterSeek -o ShatterSeek            # then: Rscript install_shatterseek_RE.R --src ShatterSeek
# one sub-directory, selected files, or a list file with one path per line
bash gh_raw_fetch.sh parklab/ShatterSeek -d R -o ShatterSeek
bash gh_raw_fetch.sh Tim-Yu/chromothripsis-shatterseek scripts/shatterseek_lib.R README.md
bash gh_raw_fetch.sh Tim-Yu/chromothripsis-shatterseek -l files.txt
# branch, tag or commit via owner/repo@ref (default main, falls back to master); a commit SHA avoids CDN cache lag
bash gh_raw_fetch.sh Tim-Yu/chromothripsis-shatterseek@d47f7c8 -o chromothripsis-shatterseek
```

## 2. Input formats

**Consensus SV bedpe** (`<PID>/<SID>/<SID>.somatic.sv.bedpe`, header optional):
```
chrom1 start1 end1 chrom2 start2 end2 sv_id pe_support strand1 strand2 svclass svmethod
```
`svclass` in `DEL DUP h2hINV t2tINV TRA` (also accepted: `INV` resolved by strand, `BND`). Strands may be
bracketed (`[+]`). VCF2BEDPE-style headers (`CHROM_A START_A ... SVTYPE`) are also recognised.

**Copy number**, one of
* Battenberg: `chr startpos endpos ... ntot|tot_n nMaj1_A nMin1_A frac1_A nMaj2_A nMin2_A frac2_A ...`
  (total CN = `ntot` if present, else `nMaj1_A+nMin1_A`, or fraction-weighted with `--subclonal weighted`)
* DRAGEN-style export: `Chromosome, Start Position, End Position, Tumour TCN, vcf_filter` (use `--filter-pass true`)
* GISTIC-style `Seg.CN` (log2 ratio, converted to absolute)

**Cohort tables** (legacy format): `SVs.txt` and `CNAs.txt` as in the table above; total CN = `nMaj1+nMin1`.

Chromosomes `chr1`/`1`/`23` are normalised to `1..22, X`; chrY and contigs are dropped (ShatterSeek supports 1–22, X).

## 3. Run

```bash
# A) per-sample layout, one sample
Rscript scripts/run_shatterseek_persample.R \
  --sv  <PID>/<SID>/<SID>.somatic.sv.bedpe \
  --cn  P_<PID>_T_<SID>_*/<SID>_copynumber_1_gapfilled_LogR_new.txt \
  --sample <SID> --patient <PID> --outdir shatterseek_out

# A) per-sample layout, batch from the matrix written by build_sv_cnv_matrix.R
#    (columns patientID sampleID SV_file_dir CN_file_dir)
Rscript scripts/run_shatterseek_persample.R --matrix sv_cnv_matrix.tsv --outdir shatterseek_out
#    --cn-col <column> to use another CN column of the matrix

# B) cohort tables
Rscript scripts/run_shatterseek_cohort.R --sv SVs.txt --cn CNAs.txt --outdir shatterseek_out_cohort \
        [--samples list.txt]

# common options
#   --genome hg38|hg19 (default hg38)   --filter-pass true|false   --subclonal clonal|weighted
#   --min-size 1   --plot true|false   --save-rds true|false
#   --thresholds k=v,k=v   e.g. --thresholds hc_min_osc2=6,p_joins=0.1  (names: default_thresholds() in shatterseek_lib.R)
```

### GEL research environment, end to end

Both GEL drivers start from a patient list (one patientID per line) and the GEL roots
(SV: `/re_gecip/cancer_sarcoma/19.ComplexSVs/19.5.SVs/19.5.7.ConsensusCalls/5_callers`,
Battenberg: `/re_gecip/cancer_sarcoma/33.CN_Sigs/33.4.BB_fix/BB_merged`), set at the top of each script.

```bash
# ShatterSeek on every sample of the patients
bash scripts/run_gel_example.sh                                   # -> shatterseek_output/shatterseek_calls.tsv

# improved legacy caller (needs the cohort tables, which are built from the same per-sample files)
bash scripts/run_gel_legacy.sh patients.txt legacy_output         # -> legacy_output/chromothripsis_improved_calls.tsv
#   also writes legacy_output/cohort_tables/{SVs.txt,CNAs.txt} and, with RUN_ORIGINAL=true (default),
#   legacy_output/chromothripsis_original_calls.tsv from the unchanged algorithm for comparison
#   LEGACY_OPTS="--minEventsSV 10" bash scripts/run_gel_legacy.sh ...   passes options to chromothripsis_improved.r
```
Either can be submitted with `bsub` (see the header of each script). The intermediate steps can be run alone:
`build_sv_cnv_matrix.R` (matrix) and `build_cohort_tables.R --matrix <tsv> --outdir <dir>` (cohort tables; strands are
derived from `svclass`, `nMaj1/nMin1/...` from the Battenberg `*_A` columns, ploidy as the length-weighted mean total CN).

### Outputs (per `--outdir`)

| File | Content |
|---|---|
| `shatterseek_calls.tsv` | one row per sample × chromosome with a call (HighConfidence / LowConfidence / Extended_ringlike) |
| `shatterseek_all_chromosomes.tsv` | all 23 chromosomes per sample with every ShatterSeek statistic and the criteria flags |
| `per_sample/<SID>_shatterseek_summary.tsv`, `per_sample/<SID>_shatterseek.rds` | per-sample table and the full ShatterSeek object |
| `plots/<SID>_chr<N>_shatterseek.pdf` | ShatterSeek SV-arc + CN plot for each called chromosome |
| `shatterseek_run_log.tsv`, `thresholds_used.txt` | status per sample, thresholds applied |

## 4. Calling rules

Published ShatterSeek criteria (package tutorial, "Criteria to call chromothripsis"):

* **HighConfidence** – (≥6 interleaved intra-chromosomal SVs, ≥7 CN segments oscillating between 2 states,
  fragment-joins test not rejected (p>0.05), and chromosome breakpoint enrichment p<0.05 or exponential
  breakpoint-distribution p<0.05) **or** (≥3 interleaved intra-chr SVs and ≥4 inter-chromosomal SVs, ≥7
  oscillating segments, fragment-joins not rejected).
* **LowConfidence** – ≥6 interleaved intra-chr SVs, 4–6 oscillating segments (2 states), fragment-joins not
  rejected, and enrichment or exponential test significant.
* **Extended_ringlike** (our addition, reported separately for manual review) – ring chromosomes / high-level
  amplicons (e.g. DFSP chr17q–chr22q rings) alternate between one baseline CN and *varying* amplified states
  (3,11,3,9,3,11,3,10,3,…), so the strict 2-state count stays small although the profile clearly oscillates.
  `max_number_oscillating_CN_segments_baseline` is the longest run in which every other segment returns to the same
  CN. Extended is called when the SV cluster including translocations has ≥10 SVs, the region oscillates
  (3-state ≥7 **or** baseline ≥7 **or** ≥10 CN segments with ≥4 translocations), and a breakpoint test is
  significant or ≥4 translocations are present.

All thresholds can be changed with `--thresholds`.

### Fixes applied to ShatterSeek 1.1 at run time (`shatterseek_lib.R`)
1. `statistical_criteria()` compares `inter$pos1` instead of `inter$pos2` when the *second* breakpoint of an
   inter-chromosomal SV lies on the candidate chromosome, so translocations stored as chrom1=partner /
   chrom2=candidate were dropped from `number_TRA` / `clusterSize_including_TRA`. The function is patched in the
   package namespace (disable with `patch = FALSE` in `run_shatterseek_sample()`).
2. The exponential-distribution test draws a random sample, so p-values were not reproducible; a seed is set per sample.

## 5. Improved legacy caller

The original PCF-density caller required ≤3 allele-specific CN states to cover 93–100 % of a region and a KS test
to reject exponential segment sizes; genuine shattered regions and multi-state ring amplicons fail both, so every
flagged region came out as "Cluster". `legacy/chromothripsis_improved.r` keeps the inputs and the 16 original
columns, makes those tests informative only, and adds oscillation statistics, an SV-breakpoint clustering test,
SV-based region flagging and linking of two-chromosome events. Calls are tiered
`Chromothripsis_HighConf` / `Chromothripsis_LowConf` / `Cluster`.

```bash
Rscript legacy/chromothripsis_improved.r --sv SVs.txt --cn CNAs.txt --out improved_calls.tsv      # --help lists ~25 thresholds
Rscript legacy/chromothripsis_original_runnable.r --sv SVs.txt --cn CNAs.txt --out original_calls.tsv
Rscript legacy/reclassify_pui_table.r <legacy_output.tsv> <reclassified.tsv>   # new rule on an existing legacy table
```
Details, thresholds and the reasoning behind each change: `legacy/README_legacy.md`.

## References
* Cortés-Ciriano I. et al. Comprehensive analysis of chromothripsis in 2,658 human cancers using whole-genome
  sequencing. *Nat Genet* 52, 331–341 (2020). https://github.com/parklab/ShatterSeek
* Korbel J.O. & Campbell P.J. Criteria for inference of chromothripsis in cancer genomes. *Cell* 152, 1226–1236 (2013).
