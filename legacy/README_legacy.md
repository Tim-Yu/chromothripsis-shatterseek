# Legacy GEL-sarcoma chromothripsis caller: algorithm, failure analysis and improved version

Files in this directory

| file | what |
|---|---|
| `chromothripsis_original_runnable.r` | same algorithm, made runnable end to end as a CLI (see section 3) |
| `chromothripsis_improved.r` | improved caller: same inputs, legacy output columns + new ones, tiered calls (section 4) |
| `reclassify_pui_table.r` | applies the new rule to an existing legacy output table (section 6) |

Commands (R >= 4.0 with `copynumber`, `MASS`, `GenomicRanges`; the mock cohort referred to below was generated locally for testing and is not part of this repository):

```bash
# legacy algorithm, unchanged logic
Rscript legacy/chromothripsis_original_runnable.r --sv mock_data/cohort/SVs.txt --cn mock_data/cohort/CNAs.txt \
        --out legacy/test_output/original_calls_mock.tsv

# improved caller (all thresholds tunable, see --help)
Rscript legacy/chromothripsis_improved.r --sv mock_data/cohort/SVs.txt --cn mock_data/cohort/CNAs.txt \
        --out legacy/test_output/improved_calls_mock.tsv
Rscript legacy/chromothripsis_improved.r --help

# apply the new final rule to an existing legacy output table
Rscript legacy/reclassify_pui_table.r Re__Pui_Boyu/All_chromothripsis_calls.txt legacy/test_output/pui_calls_reclassified.tsv
```

Input formats (identical for both scripts; this is the "Atef" cohort format):

```
SVs.txt  : samplename chr1 strand1 pos1 chr2 strand2 pos2 svclass      svclass in DEL, DUP, h2hINV, t2tINV, TRA
CNAs.txt : samplename chr start end nMaj1 nMin1 frac1 nMaj2 nMin2 frac2 SD ploidy   (Battenberg-style, one row per CN segment,
                                                                                     gap-filled so that segments tile the chromosome)
```

---

## 1. The legacy algorithm, step by step

Everything is done per sample and per chromosome (1..22, X) on the **copy-number segment table**; SVs are only used at the very end for annotation.

### Step 1 - "probes" = CN segment ends, signal = inter-breakpoint distance (`doPCF`)

For a chromosome with `n` CN segments (requires `n > minEvents = 5`, otherwise the chromosome is skipped), the script builds a
pseudo copy-number track with one probe per CN segment:

* probe position = segment `end`
* probe value = distance to the previous breakpoint: `c(end[1]-start[1], diff(end))`

and runs `copynumber::pcf(..., kmin = minEvents = 5, gamma = 40 (default), normalize = TRUE (default), assembly = "hg19" (default))`.
PCF therefore segments the chromosome into runs of probes with similar inter-breakpoint distance; each PCF segment has `n.probes` (number of CN
segments it contains), `mean` (mean inter-breakpoint distance in bp) and `arm` (p/q from the hg19 centromere table).

### Step 2 - density flag and per-arm merge (`treatPCF`, `CP`)

A PCF segment is **flagged** when

```
mean < flagDensity   with flagDensity = 150e6/50 = 3,000,000 bp      (mean inter-breakpoint distance below 3 Mb)
AND n.probes >= minEvents = 5                                        (at least 5 CN segments)
```

All flagged PCF segments of one arm are merged into a single region (from the first to the last flagged segment, including any
unflagged segments in between), so there is at most one flagged region per arm. The region reported in the table is
`start = start of the first CN segment` ... `end = end of the last CN segment` of the merged PCF segment (so it includes the telomeric
flank segment when the PCF segment starts at the first segment of the chromosome, e.g. MOCK3 `chr5p:0.0-45.0Mb`).

`density` = (probe-weighted) mean inter-breakpoint distance, `NbreakpointsCNA` = `n.probes` = **number of CN segments** in the region
(not strictly the number of breakpoints; it is breakpoints + 1 within the region).

### Step 3 - Kolmogorov-Smirnov test against an exponential (`testExp`)

`sizes = diff(end)` of the CN segments in the region (needs >= 4 sizes); an exponential is fitted by `MASS::fitdistr` and
`ks.test(sizes, pexp, rate = fit)` gives `p.KS_exp`. The idea: random breakpoints along a chromosome give exponentially distributed
gaps; a *small* p ("not exponential") was taken as evidence of non-random clustering. The same test is also run on the whole chromosome
(`testExpChromosome`) and used as an alternative flag (`flagged = tE<0.05 | flaggedDensity`) but `deriveTable()` only outputs
density-flagged regions, so the chromosome-level test has no effect on the output.

### Step 4 - copy-number state coverage (`getCNstates`, `returnCov`)

For the CN segments in the region, state = `"nMaj1-nMin1"`; the total size of each state is computed and

* `coverage.mode` = fraction of the region covered by the <= 3 states with the **most segments**
* `coverage.max`  = fraction of the region covered by the <= 3 **largest** states

`expectedCov = returnCov(NbreakpointsCNA)` is a linear function of the number of CN segments: `1.1 - 0.006 * nb`, capped at 1, i.e.

| nb CN segments | 5-16 | 17 | 22 | 29 | 50 | 100 |
|---|---|---|---|---|---|---|
| expectedCov | **1.000** | 0.998 | 0.968 | 0.926 | 0.800 | 0.500 |

The final rule requires `coverage.mode > expectedCov` **or** `coverage.max > expectedCov` (strict). Since coverage can never exceed 1,
**any region with <= 16 CN segments can never be called chromothripsis**, whatever the other statistics (this makes the `NbreakpointsCNA > 15` clause
in step 6 almost redundant). Regions with 17-30 segments need >= 93-99.8 % of their length in <= 3 allele-specific states.

### Step 5 - SV annotation (`tests`, `testclass`, `testorder`)

SVs of the sample with **at least one breakpoint inside the region** (`GenomicRanges::findOverlaps`, any chromosome for the mate) are
selected; `nbSVs` = their number. `pRandomClass` = chi-square p-value of the counts of `DEL / h2hINV / t2tINV / DUP` against
equiprobability (1/4 each; TRA are ignored; `NA` when none of the four classes is present, i.e. chisq.test errors). `pRandomOrder`
(Spearman correlation of `pos1` vs `pos2` of intrachromosomal SVs) is computed but **not used**.

### Step 6 - final rule (`isChromothripsis`)

```
Chromothripsis  <=>  (coverage.mode > expectedCov | coverage.max > expectedCov)     # <= 3 CN states cover "enough" of the region
                  &  pRandomClass > 0.01                                            # SV classes not significantly non-uniform
                  &  p.KS_exp < 0.05                                                # segment sizes "not exponential"
                  &  (NbreakpointsCNA > 15 | nbSVs > 15)                            # size of the event
otherwise         "Cluster"
```
`ifelse(NA,...)` gives `NA` when `pRandomClass` is `NA`; Pui's table shows "Cluster" for those rows, so the runnable copy maps `NA` to "Cluster".

Output columns: `samplename start end chrArm density NbreakpointsCNA p.KS_exp CNA.states CNA.sizes coverage.mode coverage.max pRandomClass pRandomOrder expectedCov nbSVs Calls`.

---

## 2. Why the DFSP (GMS) events were missed - root causes

Verified on Pui's table (`Re__Pui_Boyu/All_chromothripsis_calls.txt`) and on the mock cohort (`test_output/original_calls_mock.tsv`,
column `legacy_fail_reason` in `test_output/pui_calls_reclassified.tsv`):

1. **The KS-exponential requirement (`p.KS_exp < 0.05`) is wrong-headed and fails on real and simulated chromothripsis.** Within a shattered
   region the breakpoints *are* approximately uniformly scattered, so the inter-breakpoint distances are close to exponential and the KS
   test does *not* reject: GMS1 chr22q p = 0.051, GMS5 chr1p 0.32, GMS5 chr17q 0.61, GMS7 chr17q 0.90; MOCK3 (textbook chromothripsis,
   43 segments, coverage = 1.0) p = 0.52. Running the improved script with `--require_ks TRUE` turns every mock call into "Cluster".
   All 10 GMS regions fail this criterion except GMS9 chr6p (HLA region artefact). The test is also only informative for *random vs. clustered*
   at chromosome level, which the density flag already established.
2. **`returnCov()` demands impossible coverage.** `expectedCov = 1` for <= 16 CN segments (impossible to beat with a strict `>`), 0.93-0.998
   for 17-29 segments. DFSP ring amplicons oscillate among many high-level states (GMS1 chr22q has 15 distinct allele-specific states,
   `coverage.max = 0.61`; GMS5 chr17q 12 states, 0.61), so "<= 3 states" can never cover 93 % of them. Even GMS5 chr1p (coverage 0.98 with
   states 2-1/3-1/4-1) failed because 0.980 < 0.998. All 10 GMS regions fail this criterion.
3. **`NbreakpointsCNA > 15 | nbSVs > 15`**: a size gate that is fine for large events but, combined with (2), makes 5-16 segment events
   uncallable even with dozens of SVs (GMS1 chr5p: 9 segments, 37 SVs).
4. **`pRandomClass`** is `NA` when a region has no DEL/DUP/INV (only TRAs or no SVs), which propagates to `NA`/"Cluster"; a region made of
   translocations is penalised, and the equiprobability test on 4 classes with few SVs is noisy anyway.
5. **Region flagging is fragile.** (a) the PCF `mean` includes the distance from the arm start to the first breakpoint, so a 15 Mb telomeric
   flank pushes 17 breakpoints over 47 Mb (2.6 Mb spacing) to a mean of 3.26 Mb > 3 Mb and the region is **not flagged at all** (MOCK2 and
   MOCK4 chr1p with the legacy script; GMS5 chr1p passed with a density of 2.72 Mb only narrowly). (b) `minEvents = 5` CN segments are needed
   per PCF segment, and SV clusters whose CN changes were smoothed away by Battenberg (e.g. GMS2 chr1p 15-60 Mb, visible in the genome plot but
   absent from the table) are never considered because SVs play no role in flagging.
6. Minor: `pRandomOrder` and the chromosome-level KS test are computed but unused; `NbreakpointsCNA` is really a segment count;
   the hardcoded hg19 arm table is used on GRCh38 data; the per-sample loop was run twice (mclapply then lapply).

---

## 3. `chromothripsis_original_runnable.r` - what was changed (plumbing only)

* removed `.libPaths()`, `install.packages("EMT")`, `library(vcd)`, `library(EMT)` - neither package is used (`fitdistr` is MASS,
  `chisq.test` is base stats), so no replacement was needed;
* command line `--sv <SVs.txt> --cn <CNAs.txt> --out <out.tsv>` (and `--help`) instead of the `SV_data`/`BB_data` globals and the hardcoded
  GEL output path; no `save(allCPs)`;
* the duplicated `mclapply(mc.cores=5)` pass was removed (plain `lapply`);
* `deriveTable()` handles "nothing flagged" (writes an empty table), `NA` in the final rule is reported as "Cluster" (matches Pui's table),
  chromosome names are carried as names rather than list indices, chi-square warnings are silenced.
* Algorithm, thresholds and output columns are unchanged.

Result on the mock cohort (4 regions flagged, **0 chromothripsis**, exactly the failure mode seen on the GMS cases):

| sample | region | CN segs | p.KS_exp | coverage.mode / max | expectedCov | nbSVs | call |
|---|---|---|---|---|---|---|---|
| MOCK1 | chr17q:62.1-83.3Mb | 20 | 0.84 | 0.38 / 0.72 | 0.98 | 42 | Cluster |
| MOCK1 | chr22q:18.1-50.8Mb | 27 | 0.97 | 0.56 / 0.68 | 0.938 | 47 | Cluster |
| MOCK3 | chr5p:0.0-45.0Mb | 43 | 0.52 | 1.00 / 1.00 | 0.842 | 28 | Cluster (fails only on KS) |
| MOCK4 | chr17q:61.7-83.3Mb | 20 | 0.91 | 0.44 / 0.59 | 0.98 | 56 | Cluster |
| MOCK2 chr1p, MOCK4 chr1p, MOCK6 chr9p | - | | | | | | not flagged |

---

## 4. `chromothripsis_improved.r` - changes and reasoning

The structure of the original is kept (PCF-based region flagging -> annotation -> final rule), the input format and the 16 legacy
columns are unchanged (`Calls` now holds the tiered call; the legacy rule is still evaluated into `Calls_legacyRule`). All thresholds are
command-line options (`--help` lists them with defaults).

### 4.1 Region flagging

| change | reasoning |
|---|---|
| flank-insensitive density: a PCF segment is flagged when its legacy `mean` **or** its span density `(end.pos-start.pos)/(n.probes-1)` is `< flagDensity` (3 Mb) | the legacy mean includes the gap from the arm start to the first breakpoint; a long telomeric flank hid MOCK2/MOCK4 chr1p (root cause 5a) |
| region coordinates = hull of the CN breakpoints inside the PCF segment (no telomeric flank) | cleaner coordinates for the SV clustering test; otherwise identical to legacy |
| **SV-based flagging** (`--flag_from_sv TRUE`, `--minEventsSV 8`): the same PCF/density procedure is run on the sorted SV breakpoint positions of the chromosome; SV-derived hulls are trimmed of isolated terminal breakpoints (`--trim_gap 6 Mb`); CN- and SV-derived regions on the same arm are merged when they overlap (`flagSource` = CN, SV or CN+SV) | SV clusters without matching CN segmentation (GMS2 chr1p) are now evaluated; they end up as "Cluster" unless CN oscillation is present |
| `--assembly hg38` for arm assignment (legacy: hg19) | GEL data are GRCh38 |

### 4.2 Annotation - legacy statistics kept, new ones added

Legacy (unchanged definitions): `density, NbreakpointsCNA (= CN segments overlapping the region), p.KS_exp, CNA.states, CNA.sizes, coverage.mode, coverage.max, pRandomClass, pRandomOrder, expectedCov, nbSVs`.

New columns:

| column | definition |
|---|---|
| `flagSource`, `region_size_Mb`, `span_density`, `nCNbp` | how the region was found; size; size / internal breakpoints; number of CN breakpoints (= segments - 1) |
| `CN_total`, `nCNstates`, `minCN`, `maxCN`, `ploidy`, `highLevelAmp` | total CN (`nMaj1+nMin1`) per segment; distinct states; `maxCN >= ploidy + 5` |
| `osc2` | ShatterSeek-like: longest run of adjacent segments alternating between exactly **2** total-CN states (`cn[i]==cn[i-2] != cn[i-1]`), in segments |
| `osc3` | longest run of adjacent-differing segments using at most **3** distinct total-CN states |
| `fracDirChange` | fraction of up/down direction reversals along the CN profile (zigzag). ~2/3 for a random walk, ~0 for a staircase (BFB, stepwise gains), ~1 for chromothripsis / ring amplicons that keep jumping between unrelated levels (GMS1 chr22q: 0.82 with 15 distinct states) |
| `CN_pattern` | `oscillating_2states` / `oscillating_2-3states` / `multistate_highlevel_amplicon` / `multistate_zigzag` / `other` |
| `nSVbp_region`, `nSV_intra`, `nSV_partialChr`, `nSV_inter`, `partnerChr`, `frac_inter_topPartner` | SV breakpoints in the region; SVs with both ends inside; one end inside + mate elsewhere on the same chromosome; inter-chromosomal; top-3 partner chromosomes with counts; fraction of inter-chromosomal SVs going to the top partner |
| `nSVbp_chr`, `p.SVcluster`, `p.SVcluster_genome` | SV breakpoints on the chromosome; one-sided binomial test of `nSVbp_region` vs. expected under a uniform distribution of the chromosome's (genome's) SV breakpoints (`p = region length / chromosome (genome) length`) |
| `fracCNbp_withSV` | fraction of CN breakpoints with an SV breakpoint within `--bp_tol 20 kb` (QC of CN/SV concordance) |
| `linkedTo`, `nSV_linked` | other flagged regions of the same sample joined by `>= --min_inter_link 5` inter-chromosomal SVs (DFSP chr17q<->chr22q rings, GMS5 chr1p<->chr17q) |
| `Calls_legacyRule`, `evidence` | the legacy rule on the same region; human-readable summary of the evidence that produced `Calls` |

### 4.3 Final rule (tiered)

```
oscStrong   = osc2 >= 7  |  osc3 >= 9  |  (fracDirChange >= 0.7 & NbreakpointsCNA >= 15)
oscWeak     = osc2 >= 4  |  osc3 >= 6  |  (fracDirChange >= 0.6 & NbreakpointsCNA >= 8)
clustStrong = nbSVs >= 15 & p.SVcluster <= 0.01
clustWeak   = nbSVs >= 8  & p.SVcluster <= 0.05
linked      = >= 5 inter-chromosomal SVs join the region to another flagged region of the sample that itself has oscWeak & clustWeak

Chromothripsis_HighConf = NbreakpointsCNA >= 10 & clustStrong & (oscStrong | (oscWeak & linked))
Chromothripsis_LowConf  = NbreakpointsCNA >= 7  & clustWeak   & oscWeak            (and not HighConf)
Cluster                 = everything else that was density-flagged
```

Optional legacy gates: `--require_ks TRUE` (p.KS_exp < 0.05) and `--require_random_class TRUE` (pRandomClass > 0.01); both default to
FALSE and the p-values are reported as information only.

Reasoning behind the replacements:

* **KS test -> informative only** (root cause 1). Clustering is instead established by the density flag plus an explicit SV-breakpoint
  clustering test (`p.SVcluster`), which is what the KS test was meant to capture.
* **`returnCov` coverage -> oscillation statistics** (root cause 2). "Few CN states covering the region" was a proxy for the
  chromothripsis hallmark *oscillation between a small number of states*; `osc2`/`osc3` measure that directly (as ShatterSeek: >= 7
  oscillating segments for high confidence, 4-6 for low confidence), and the zigzag fraction captures high-level multi-state amplicons
  (DFSP rings) that oscillate but never settle on <= 3 states. Coverage columns are still reported.
* **Hard size gate -> graded gates** (root cause 3): 10 CN segments & 15 SVs for HighConf, 7 & 8 for LowConf. Small sparse clusters
  (GMS6/GMS8 chr9p, MOCK6) stay "Cluster".
* **SV support is mandatory** for a chromothripsis call (>= 8 SVs with a breakpoint in the region). This is what stops segmentation
  artefacts with beautiful oscillation but no SVs (GMS9 chr6p HLA region: osc2 = 5, osc3 = 8, 0 SVs) from being called.
* **pRandomClass** no longer kills translocation-dominated regions (root cause 4); inter-chromosomal SVs are counted (`nSV_inter`) and used
  positively through region linking, which recognises two-chromosome events (DFSP chr17q/chr22q, GMS5 chr1p/chr17q).
* **Tiered output**: HighConf requires strong oscillation or a linked partner region; LowConf flags events worth looking at in the genome plot.

---

## 5. Mock cohort: truth vs original vs improved (`test_output/mock_comparison.tsv`)

| sample | planted event | expected | original (region / call) | improved (region / call) | OK |
|---|---|---|---|---|---|
| MOCK1 | chr22:15.9-50.7Mb DFSP_ring_amplicon | chromothripsis-like | chr22q:18.1-50.8Mb / Cluster | chr22q:15.9-50.7Mb / Chromothripsis_HighConf | yes |
| MOCK1 | chr17:56.2-83.2Mb DFSP_ring_amplicon | chromothripsis-like | chr17q:62.1-83.3Mb / Cluster | chr17q:56.2-83.2Mb / Chromothripsis_HighConf | yes |
| MOCK2 | chr1:15.0-62.0Mb clustered_rearrangement_oscillating_gain | chromothripsis-like | not flagged / - | chr1p:15.0-62.0Mb / Chromothripsis_HighConf | yes |
| MOCK3 | chr5:1.0-45.0Mb canonical_chromothripsis | high_confidence | chr5p:0.0-45.0Mb / Cluster | chr5p:1.0-45.0Mb / Chromothripsis_HighConf | yes |
| MOCK4 | chr1:16.6-62.9Mb two_chromosome_amplicon | chromothripsis-like | not flagged / - | chr1p:16.6-62.9Mb / Chromothripsis_HighConf | yes |
| MOCK4 | chr17:56.2-83.2Mb two_chromosome_amplicon | chromothripsis-like | chr17q:61.7-83.3Mb / Cluster | chr17q:56.2-83.2Mb / Chromothripsis_HighConf | yes |
| MOCK5 | none (negative control) | none | nothing flagged / - | nothing flagged / - | yes |
| MOCK6 | chr9:39.0-43.2Mb small_cluster | cluster_only | not flagged / - | chr9p:39.0-42.8Mb / Cluster | yes |

No region outside the planted events was flagged. Key evidence per detected region (from `improved_calls_mock.tsv`):

| region | CN segs | nbSVs (intra / inter) | osc2 / osc3 / zigzag | p.SVcluster | linked | pattern |
|---|---|---|---|---|---|---|
| MOCK1 chr22q | 32 | 52 (22 / 30 to chr17) | 5 / 7 / 0.93 | 7e-13 | chr17q:30 | oscillating_2-3states |
| MOCK1 chr17q | 25 | 49 (19 / 30 to chr22) | 3 / 5 / 0.91 | 6e-34 | chr22q:30 | multistate_highlevel_amplicon |
| MOCK2 chr1p | 19 | 69 (61 / 8) | 3 / 5 / 1.00 | 8e-93 | - | multistate_zigzag |
| MOCK3 chr5p | 43 | 28 (28 / 0) | 32 / 32 / 1.00 | 3e-35 | - | oscillating_2states |
| MOCK4 chr1p | 19 | 76 (25 / 52 to chr17) | 3 / 5 / 1.00 | 2e-69 | chr17q:50 | multistate_zigzag |
| MOCK4 chr17q | 25 | 75 (25 / 50 to chr1) | 0 / 5 / 0.83 | 1e-49 | chr1p:50 | multistate_highlevel_amplicon |
| MOCK6 chr9p | 6 | 3 | 4 / 4 / 1.00 | 1e-8 | - | Cluster (too few CN segments and SVs) |

Sensitivity checks: with `--flag_from_sv FALSE` the same 7 regions are flagged (the span-density fix alone recovers MOCK2/MOCK4 chr1p);
with `--require_ks TRUE` all 7 become "Cluster" (the KS criterion is the single most damaging legacy requirement).

Caveat: the mock SV breakpoints coincide exactly with CN breakpoints and intra/inter SVs are drawn randomly; real data will have noisier
`fracCNbp_withSV` and `osc` values. Thresholds are deliberately conservative on the SV side (>= 15 SVs for HighConf) to compensate.

---

## 6. Pui's 10 GMS regions under the new rule (`test_output/pui_calls_reclassified.tsv`)

The legacy table carries `CNA.states`, so total CN per segment, `osc2`, `osc3` and `fracDirChange` can be recomputed exactly, together with
`NbreakpointsCNA` and `nbSVs`. **Not computable from the table** (needs raw SV coordinates): `p.SVcluster` (assumed passed because every
row is a density-flagged region - the only risk is a region with many SVs that are *not* clustered within it), `nSV_intra / nSV_inter /
partnerChr`, `linkedTo` (so the "linked" upgrade LowConf -> HighConf is not applied), `fracCNbp_withSV`, `highLevelAmp` (no ploidy column).
Rerunning `chromothripsis_improved.r` on Pui's `SVs.txt`/`CNAs.txt` would fill these in.

| GMS | region | CNbp | SVs | CN_total | osc2 | osc3 | zigzag | legacy call | legacy failed on | new call |
|---|---|---|---|---|---|---|---|---|---|---|
| GMS2 | chr1p:107.5-122.0Mb | 5 | 7 | 4;5;4;5;3 | 4 | 5 | 1 | Cluster | coverage; KS; size | **Cluster** |
| GMS1 | chr5p:0.0-14.5Mb | 9 | 37 | 8;12;4;10;4;6;6;9;7 | 3 | 4 | 0.83 | Cluster | coverage; KS | **Chromothripsis_LowConf** |
| GMS1 | chr22q:15.9-50.8Mb | 29 | 62 | 5;1;4;2;3;16;5;10;5;3;1;5;3;11;3;9;3;11;3;10;3;4;1;8;9;7;10;7;0 | 3 | 7 | 0.82 | Cluster | coverage; KS | **Chromothripsis_HighConf** |
| GMS8 | chr9p:39.0-42.5Mb | 5 | 0 | 3;4;5;4;3 | 3 | 5 | 0.33 | Cluster | coverage; class NA; KS; size | **Cluster** |
| GMS5 | chr1p:16.6-62.9Mb | 17 | 75 | 15;4;3;3;7;3;5;3;3;4;6;4;4;5;4;3;3 | 3 | 5 | 0.73 | Cluster | coverage; KS | **Chromothripsis_HighConf** |
| GMS5 | chr17q:56.2-83.2Mb | 22 | 77 | 5;4;8;9;11;9;8;10;8;7;6;8;10;6;8;10;6;7;9;7;9;7 | 5 | 7 | 0.60 | Cluster | coverage; KS | **Chromothripsis_LowConf** (would become HighConf through the chr1p link if the raw SVs confirm >= 5 chr1p<->chr17q joins, as the genome plot suggests) |
| GMS5 | chr21p:5.2-10.8Mb | 6 | 1 | 15;12;7;13;25;7 | 0 | 4 | 0.5 | Cluster | coverage; KS; size | **Cluster** |
| GMS9 | chr6p:29.9-32.7Mb | 8 | 0 | 1;0;2;0;1;0;1;0 | 5 | 8 | 1 | Cluster | coverage; class NA; size | **Cluster** (oscillates, but 0 SVs: HLA-region segmentation artefact) |
| GMS6 | chr9p:39.0-43.2Mb | 5 | 0 | 5;3;6;0;2 | 0 | 3 | 1 | Cluster | coverage; class NA; KS; size | **Cluster** |
| GMS7 | chr17q:52.8-59.5Mb | 12 | 25 | 3;5;3;5;3;5;3;5;3;5;3;4 | 11 | 12 | 1 | Cluster | coverage; KS | **Chromothripsis_HighConf** (11 segments alternating 3<->5 with 25 SVs; COL1A1 region joined to chr22 in the genome plot) |

Summary: 5 of the 10 regions flip - GMS1 chr22q, GMS5 chr1p and GMS7 chr17q to **HighConf**; GMS1 chr5p and GMS5 chr17q to **LowConf**.
The five regions with <= 8 CN segments and <= 7 SVs stay "Cluster". GMS2's chr1p 15-60 Mb SV cluster is absent from Pui's table (never
flagged by the CN-only legacy procedure) and can only be assessed by rerunning the improved script, whose SV-based flagging targets exactly
that case.
