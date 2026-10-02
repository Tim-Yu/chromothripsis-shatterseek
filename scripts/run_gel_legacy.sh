#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_gel_legacy.sh - run the (improved) legacy chromothripsis caller in the GEL
# research environment on the per-sample consensus SV / Battenberg files.
#
#   step 1  patient list -> matrix of SV bedpe + Battenberg CN paths  (build_sv_cnv_matrix.R)
#   step 2  matrix       -> cohort tables SVs.txt + CNAs.txt            (build_cohort_tables.R)
#   step 3  cohort tables -> legacy/chromothripsis_improved.r            (+ the unchanged original for comparison)
#
# Usage:
#   bash run_gel_legacy.sh [patients.txt] [outdir]
#   bsub -q medium -P re_gecip_cancer_sarcoma -n 2 -R "rusage[mem=8000]" \
#        -o legacy.%J.log -e legacy.%J.err bash run_gel_legacy.sh patients.txt legacy_output
#
# Edit the path variables below if the GEL layout changes. Extra options for the
# improved caller can be passed through LEGACY_OPTS, e.g.
#   LEGACY_OPTS="--minEventsSV 10 --assembly hg38" bash run_gel_legacy.sh
# ---------------------------------------------------------------------------
set -euo pipefail

module load R/4.2.1 2>/dev/null || true

PATIENTS_TXT="${1:-patients.txt}"                                 # one patientID per line
OUTDIR="${2:-legacy_output}"
SV_ROOT="/re_gecip/cancer_sarcoma/19.ComplexSVs/19.5.SVs/19.5.7.ConsensusCalls/5_callers"
CN_ROOT="/re_gecip/cancer_sarcoma/33.CN_Sigs/33.4.BB_fix/BB_merged"   # holds P_<PID>_T_<SID>_* folders
RUN_ORIGINAL="${RUN_ORIGINAL:-true}"                              # also run the unchanged legacy algorithm
LEGACY_OPTS="${LEGACY_OPTS:-}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LEGACY_DIR="$(cd "${HERE}/../legacy" && pwd)"
mkdir -p "${OUTDIR}"

echo "== step 1: matrix"
Rscript "${HERE}/build_sv_cnv_matrix.R" \
  --patients "${PATIENTS_TXT}" --sv-root "${SV_ROOT}" --cn-root "${CN_ROOT}" \
  --out "${OUTDIR}/sv_cnv_matrix.tsv"
#  defaults: --sv-file-pattern '{sampleID}\.somatic\.sv\.bedpe'
#            --cn-dir-pattern  'P_{patientID}_T_{sampleID}_.*'
#            --cn-file-pattern '{sampleID}_copynumber_1_gapfilled_LogR_new\.txt'

echo "== step 2: cohort tables"
Rscript "${HERE}/build_cohort_tables.R" --matrix "${OUTDIR}/sv_cnv_matrix.tsv" --outdir "${OUTDIR}/cohort_tables"

echo "== step 3: improved legacy caller"
# shellcheck disable=SC2086
Rscript "${LEGACY_DIR}/chromothripsis_improved.r" \
  --sv "${OUTDIR}/cohort_tables/SVs.txt" --cn "${OUTDIR}/cohort_tables/CNAs.txt" \
  --out "${OUTDIR}/chromothripsis_improved_calls.tsv" ${LEGACY_OPTS}

if [ "${RUN_ORIGINAL}" = "true" ]; then
  echo "== step 3b: original legacy algorithm (for comparison)"
  Rscript "${LEGACY_DIR}/chromothripsis_original_runnable.r" \
    --sv "${OUTDIR}/cohort_tables/SVs.txt" --cn "${OUTDIR}/cohort_tables/CNAs.txt" \
    --out "${OUTDIR}/chromothripsis_original_calls.tsv"
fi

echo "Calls: ${OUTDIR}/chromothripsis_improved_calls.tsv"
# The same cohort tables can also be fed to ShatterSeek:
#   Rscript "${HERE}/run_shatterseek_cohort.R" --sv "${OUTDIR}/cohort_tables/SVs.txt" --cn "${OUTDIR}/cohort_tables/CNAs.txt" --outdir "${OUTDIR}/shatterseek"
