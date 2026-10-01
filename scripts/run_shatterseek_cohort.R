#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# run_shatterseek_cohort.R
#
# ShatterSeek chromothripsis detection for the cohort-table layout used by the
# legacy chromothripsis.r (Dan / Atef):
#
#   SVs.txt  : samplename chr1 strand1 pos1 chr2 strand2 pos2 svclass
#              (svclass in DEL, DUP, h2hINV, t2tINV, TRA; strands + / -)
#   CNAs.txt : samplename chr start end nMaj1 nMin1 frac1 nMaj2 nMin2 frac2 SD ploidy
#              (Battenberg clonal/subclonal states collapsed per sample)
#
# In the GEL RE these were ../INPUT/SVs.txt and ../INPUT/CNAs.txt relative to
# ~/re_gecip/cancer_sarcoma/3.landscape/3.4.SVs/3.4.4.Chromothripsis_calls/
#
# Usage:
#   Rscript run_shatterseek_cohort.R --sv <SVs.txt> --cn <CNAs.txt> --outdir <dir>
#          [--samples <file with one samplename per line>] [--genome hg38] [--subclonal clonal|weighted]
#          [--min-size 1] [--plot true] [--save-rds true] [--thresholds k=v,...]
#
# Outputs: same as run_shatterseek_persample.R (per_sample/, plots/, shatterseek_calls.tsv, ...)
# ---------------------------------------------------------------------------

script_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})()
source(file.path(script_dir, "shatterseek_lib.R"))

usage <- function() cat(
  "Usage:\n  Rscript run_shatterseek_cohort.R --sv <SVs.txt> --cn <CNAs.txt> --outdir <dir>\n",
  "       [--samples <list>] [--genome hg38] [--subclonal clonal] [--min-size 1] [--plot true] [--save-rds true] [--thresholds k=v,...]\n", sep = "")

main <- function() {
  args <- parse_cli()
  if (isTRUE(args$help) || length(args) == 0) { usage(); return(invisible()) }
  sv_path <- args$sv %||% stop("Missing --sv")
  cn_path <- args$cn %||% stop("Missing --cn")
  outdir <- args$outdir %||% stop("Missing --outdir")
  genome <- args$genome %||% "hg38"
  subclonal <- args$subclonal %||% "clonal"
  min_size <- as.integer(args$`min-size` %||% "1")
  do_plot <- parse_bool(args$plot, TRUE)
  save_rds <- parse_bool(args$`save-rds`, TRUE)
  th <- parse_thresholds(args$thresholds)

  svs <- read_cohort_sv(sv_path)
  cnas <- read_cohort_cn(cn_path, subclonal = subclonal)
  samples <- unique(as.character(cnas$samplename))
  if (!is.null(args$samples)) {
    want <- trimws(readLines(args$samples)); want <- want[nzchar(want)]
    missing <- setdiff(want, samples); if (length(missing)) warning("Not in CNAs.txt: ", paste(missing, collapse = ", "))
    samples <- intersect(samples, want)
  }
  no_sv <- setdiff(samples, unique(as.character(svs$samplename)))
  if (length(no_sv)) { message("No SVs for ", length(no_sv), " sample(s); skipped: ", paste(head(no_sv, 10), collapse = ", ")); samples <- setdiff(samples, no_sv) }

  sample_df <- data.frame(patientID = NA, sampleID = samples, stringsAsFactors = FALSE)
  load_fun <- function(i) {
    sid <- samples[i]
    s <- svs[as.character(svs$samplename) == sid, ]
    c <- cnas[as.character(cnas$samplename) == sid, ]
    list(sv = data.frame(chrom1 = s$chr1, pos1 = s$pos1, chrom2 = s$chr2, pos2 = s$pos2, svclass = s$svclass,
                         strand1 = s$strand1, strand2 = s$strand2, stringsAsFactors = FALSE),
         cn = data.frame(chrom = c$chr, start = c$start, end = c$end, total_cn = c$total_cn, stringsAsFactors = FALSE))
  }
  message("Genome: ", genome, " | samples: ", length(samples), " | outdir: ", outdir)
  run_shatterseek_batch(sample_df, load_fun, outdir, genome = genome, th = th, min_size = min_size, do_plot = do_plot, save_rds = save_rds)
}

main()
