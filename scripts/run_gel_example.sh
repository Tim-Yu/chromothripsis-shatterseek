#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_gel_example.sh - end-to-end example for the GEL research environment.
#
# Step 1 builds the patientID/sampleID -> SV bedpe / Battenberg CN matrix with
#        the group's existing build_sv_cnv_matrix.R (from the circos_plot folder).
# Step 2 runs ShatterSeek on every row of that matrix.
#
# Edit the four path variables below for your environment, then:
#   bash run_gel_example.sh                      # run interactively
#   bsub -q medium -P re_gecip_cancer_sarcoma -n 4 -R "rusage[mem=8000]" \
#        -o shatterseek.%J.log -e shatterseek.%J.err bash run_gel_example.sh   # or via LSF
#
# ShatterSeek is not on CRAN. Inside the RE (no GitHub access) install it from
# the tarball shipped in ../deps:
#   R CMD INSTALL ../deps/ShatterSeek-master.tar.gz
# or  Rscript -e 'install.packages("../deps/ShatterSeek-master.tar.gz", repos=NULL, type="source")'
# Dependencies (BiocGenerics, graph, S4Vectors, GenomicRanges, IRanges, MASS,
# ggplot2, gridExtra, foreach) are in the standard R modules.
# ---------------------------------------------------------------------------
set -euo pipefail

module load R/4.2.1 2>/dev/null || true

# --- paths -----------------------------------------------------------------
PATIENTS_TXT="patients.txt"                                       # one patientID per line
SV_ROOT="/re_gecip/cancer_sarcoma/19.ComplexSVs/19.5.SVs/19.5.7.ConsensusCalls/5_callers"
CN_ROOT="/re_gecip/cancer_sarcoma/<battenberg_output_root>"       # holds P_<PID>_T_<SID>_* folders
CIRCOS_SCRIPTS="/re_gecip/cancer_sarcoma/32.TERT_project/Pui_BEDfiles/circos_plot"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="shatterseek_output"
GENOME="hg38"

mkdir -p "${OUTDIR}"

# --- step 1: matrix ---------------------------------------------------------
Rscript "${CIRCOS_SCRIPTS}/build_sv_cnv_matrix.R" \
  --patients "${PATIENTS_TXT}" \
  --sv-root "${SV_ROOT}" \
  --cn-root "${CN_ROOT}" \
  --out "${OUTDIR}/sv_cnv_matrix.tsv"
#  defaults: --sv-file-pattern '{sampleID}\.somatic\.sv\.bedpe'
#            --cn-dir-pattern  'P_{patientID}_T_{sampleID}_.*'
#            --cn-file-pattern '{sampleID}_copynumber_1_gapfilled_LogR_new\.txt'

# --- step 2: ShatterSeek ----------------------------------------------------
Rscript "${HERE}/run_shatterseek_persample.R" \
  --matrix "${OUTDIR}/sv_cnv_matrix.tsv" \
  --outdir "${OUTDIR}" \
  --genome "${GENOME}" \
  --subclonal clonal \
  --plot true

echo "Calls: ${OUTDIR}/shatterseek_calls.tsv"

# --- alternative: the legacy cohort tables -----------------------------------
# Rscript "${HERE}/run_shatterseek_cohort.R" \
#   --sv ~/re_gecip/cancer_sarcoma/3.landscape/3.4.SVs/INPUT/SVs.txt \
#   --cn ~/re_gecip/cancer_sarcoma/3.landscape/3.4.SVs/INPUT/CNAs.txt \
#   --outdir "${OUTDIR}_cohort" --genome hg38
