#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# run_shatterseek_persample.R
#
# ShatterSeek chromothripsis detection for the per-sample file layout used by
# the circos pipeline (circos_plot_new.R / build_sv_cnv_matrix.R):
#
#   SV : <sv_root>/<patientID>/<sampleID>/<sampleID>.somatic.sv.bedpe
#        (consensus-caller bedpe; columns chrom1 start1 end1 chrom2 start2 end2
#         sv_id pe_support strand1 strand2 svclass svmethod)
#        default sv_root in the GEL RE:
#        /re_gecip/cancer_sarcoma/19.ComplexSVs/19.5.SVs/19.5.7.ConsensusCalls/5_callers/
#   CN : Battenberg  <cn_root>/P_<patientID>_T_<sampleID>_*/<sampleID>_copynumber_1_gapfilled_LogR_new.txt
#        (chr startpos endpos ... ntot|tot_n nMaj1_A nMin1_A frac1_A ...)
#        or a DRAGEN-style export (Chromosome, Start Position, End Position, Tumour TCN, vcf_filter)
#
# Two ways to run:
#   (1) one sample:
#       Rscript run_shatterseek_persample.R --sv <bedpe> --cn <cn_file> --sample <SID> [--patient <PID>] --outdir <dir>
#   (2) a batch from the matrix written by build_sv_cnv_matrix.R
#       (columns patientID sampleID SV_file_dir CN_file_dir):
#       Rscript run_shatterseek_persample.R --matrix <matrix.tsv> --outdir <dir> [--cn-col CN_file_dir]
#
# Options:
#   --genome hg38|hg19        (default hg38)
#   --filter-pass true|false  keep only CN rows with vcf_filter == PASS (DRAGEN files; default false)
#   --subclonal clonal|weighted  total CN from Battenberg: clonal state (nMaj1+nMin1) or fraction-weighted
#   --min-size <int>          minimum interleaved-cluster size passed to shatterseek() (default 1)
#   --plot true|false         draw ShatterSeek SV/CN plots for called chromosomes (default true)
#   --save-rds true|false     save the full ShatterSeek object per sample (default true)
#   --thresholds k=v,k=v      override calling thresholds (see default_thresholds() in shatterseek_lib.R)
#
# Outputs (under --outdir):
#   per_sample/<SID>_shatterseek_summary.tsv   all 23 chromosomes with statistics + call
#   per_sample/<SID>_shatterseek.rds           ShatterSeek object (for plot_chromothripsis etc.)
#   plots/<SID>_chr<N>_shatterseek.pdf         for called chromosomes
#   shatterseek_all_chromosomes.tsv            concatenated per-chromosome table
#   shatterseek_calls.tsv                      only chromosomes with a call
#   shatterseek_run_log.tsv, thresholds_used.txt
# ---------------------------------------------------------------------------

script_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})()
source(file.path(script_dir, "shatterseek_lib.R"))

usage <- function() cat(
  "Usage:\n",
  "  Rscript run_shatterseek_persample.R --sv <bedpe> --cn <cn_file> --sample <SID> [--patient <PID>] --outdir <dir> [options]\n",
  "  Rscript run_shatterseek_persample.R --matrix <sv_cnv_matrix.tsv> --outdir <dir> [--cn-col CN_file_dir] [options]\n",
  "Options: --genome hg38 --filter-pass false --subclonal clonal --min-size 1 --plot true --save-rds true --thresholds k=v,...\n",
  sep = "")

main <- function() {
  args <- parse_cli()
  if (isTRUE(args$help) || length(args) == 0) { usage(); return(invisible()) }
  outdir <- args$outdir %||% stop("Missing --outdir")
  genome <- args$genome %||% "hg38"
  filter_pass <- parse_bool(args$`filter-pass`, FALSE)
  subclonal <- args$subclonal %||% "clonal"
  min_size <- as.integer(args$`min-size` %||% "1")
  do_plot <- parse_bool(args$plot, TRUE)
  save_rds <- parse_bool(args$`save-rds`, TRUE)
  th <- parse_thresholds(args$thresholds)

  if (!is.null(args$matrix)) {
    m <- read.delim(args$matrix, check.names = FALSE, stringsAsFactors = FALSE)
    cn_col <- args$`cn-col` %||% "CN_file_dir"
    need <- c("patientID", "sampleID", "SV_file_dir", cn_col)
    miss <- setdiff(need, colnames(m)); if (length(miss)) stop("matrix missing columns: ", paste(miss, collapse = ", "))
    ok <- !is.na(m$SV_file_dir) & nzchar(m$SV_file_dir) & !is.na(m[[cn_col]]) & nzchar(m[[cn_col]])
    if (any(!ok)) message("Skipping ", sum(!ok), " row(s) without SV or CN path")
    m <- m[ok, , drop = FALSE]
    samples <- data.frame(patientID = m$patientID, sampleID = m$sampleID, sv = m$SV_file_dir, cn = m[[cn_col]], stringsAsFactors = FALSE)
  } else {
    sv <- args$sv %||% stop("Missing --sv (or --matrix)")
    cn <- args$cn %||% stop("Missing --cn (or --matrix)")
    sid <- args$sample %||% sub("\\.somatic.*$|\\.bedpe$|\\.tsv$", "", basename(sv))
    samples <- data.frame(patientID = args$patient %||% NA, sampleID = sid, sv = sv, cn = cn, stringsAsFactors = FALSE)
  }

  load_fun <- function(i) {
    list(sv = read_sv_bedpe(samples$sv[i]),
         cn = read_cn_file(samples$cn[i], filter_pass = filter_pass, subclonal = subclonal))
  }
  message("Genome: ", genome, " | samples: ", nrow(samples), " | outdir: ", outdir)
  run_shatterseek_batch(samples, load_fun, outdir, genome = genome, th = th, min_size = min_size, do_plot = do_plot, save_rds = save_rds)
}

main()
