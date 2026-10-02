#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# build_cohort_tables.R
#
# Converts the per-sample GEL files (consensus SV bedpe + Battenberg copy
# number, listed in the matrix written by build_sv_cnv_matrix.R) into the two
# cohort tables expected by the legacy callers (legacy/chromothripsis_*.r) and
# by run_shatterseek_cohort.R:
#
#   SVs.txt  : samplename chr1 strand1 pos1 chr2 strand2 pos2 svclass
#   CNAs.txt : samplename chr start end nMaj1 nMin1 frac1 nMaj2 nMin2 frac2 SD ploidy
#
# Usage:
#   Rscript build_cohort_tables.R --matrix <sv_cnv_matrix.tsv> --outdir <dir>
#          [--cn-col CN_file_dir] [--filter-pass true|false] [--sample-col sampleID]
#
# Notes
#   * samplename = sampleID column of the matrix (e.g. LP3001668-DNA_A02)
#   * strands of intra-chromosomal SVs are derived from svclass
#     (DEL +/-, DUP -/+, h2hINV +/+, t2tINV -/-); TRA strands are taken from the bedpe
#   * nMaj1/nMin1/frac1/nMaj2/nMin2/frac2/SD come from the Battenberg columns
#     nMaj1_A ... SDfrac_A; if the file only carries a total CN (ntot / Tumour TCN)
#     nMaj1 = ceiling(cn/2), nMin1 = floor(cn/2), frac1 = 1
#   * ploidy = length-weighted mean total CN over chromosomes 1-22 (if no ploidy column)
#   * chromosome names are written without "chr" (1..22, X), as in the legacy tables
# ---------------------------------------------------------------------------

script_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})()
source(file.path(script_dir, "shatterseek_lib.R"))

main <- function() {
  args <- parse_cli()
  if (isTRUE(args$help) || length(args) == 0) {
    cat("Usage: Rscript build_cohort_tables.R --matrix <sv_cnv_matrix.tsv> --outdir <dir> [--cn-col CN_file_dir] [--filter-pass false] [--sample-col sampleID]\n")
    return(invisible())
  }
  matrix_file <- args$matrix %||% stop("Missing --matrix")
  outdir <- args$outdir %||% stop("Missing --outdir")
  cn_col <- args$`cn-col` %||% "CN_file_dir"
  sample_col <- args$`sample-col` %||% "sampleID"
  filter_pass <- parse_bool(args$`filter-pass`, FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  m <- read.delim(matrix_file, check.names = FALSE, stringsAsFactors = FALSE)
  need <- c(sample_col, "SV_file_dir", cn_col)
  miss <- setdiff(need, colnames(m)); if (length(miss)) stop("matrix missing columns: ", paste(miss, collapse = ", "))
  ok <- !is.na(m$SV_file_dir) & nzchar(m$SV_file_dir) & !is.na(m[[cn_col]]) & nzchar(m[[cn_col]])
  if (any(!ok)) message("Skipping ", sum(!ok), " row(s) without SV or CN path")
  m <- m[ok, , drop = FALSE]

  sv_rows <- list(); cn_rows <- list(); log <- list()
  for (i in seq_len(nrow(m))) {
    sid <- m[[sample_col]][i]
    res <- try({
      # --- SVs --------------------------------------------------------------
      sv <- read_sv_bedpe(m$SV_file_dir[i])
      sv$chrom1 <- normalize_chr_bare(sv$chrom1); sv$chrom2 <- normalize_chr_bare(sv$chrom2)
      sv$strand1 <- clean_strand(sv$strand1); sv$strand2 <- clean_strand(sv$strand2)
      same <- sv$chrom1 == sv$chrom2
      sv$svclass <- canonical_svtype(sv$svclass, sv$strand1, sv$strand2, same)
      for (k in names(CLASS_STRANDS)) {
        idx <- which(sv$svclass == k)
        sv$strand1[idx] <- CLASS_STRANDS[[k]][1]; sv$strand2[idx] <- CLASS_STRANDS[[k]][2]
      }
      sv$strand1[is.na(sv$strand1)] <- "+"; sv$strand2[is.na(sv$strand2)] <- "-"
      keep <- sv$chrom1 %in% SHATTERSEEK_CHROMS & sv$chrom2 %in% SHATTERSEEK_CHROMS & !is.na(sv$svclass) & !is.na(sv$pos1) & !is.na(sv$pos2)
      sv <- sv[keep, ]
      swap <- sv$chrom1 == sv$chrom2 & sv$pos1 > sv$pos2
      if (any(swap)) { tmp <- sv$pos1[swap]; sv$pos1[swap] <- sv$pos2[swap]; sv$pos2[swap] <- tmp }
      sv_rows[[sid]] <- data.frame(samplename = sid, chr1 = sv$chrom1, strand1 = sv$strand1, pos1 = sv$pos1,
                                   chr2 = sv$chrom2, strand2 = sv$strand2, pos2 = sv$pos2, svclass = sv$svclass,
                                   stringsAsFactors = FALSE)

      # --- CN ---------------------------------------------------------------
      raw <- read.delim(m[[cn_col]][i], check.names = FALSE, stringsAsFactors = FALSE)
      raw <- resolve_cn_cols(raw)
      if (filter_pass && "vcf_filter" %in% colnames(raw)) raw <- raw[tolower(raw$vcf_filter) == "pass", , drop = FALSE]
      raw$chrom <- normalize_chr_bare(raw$chrom)
      raw <- raw[raw$chrom %in% SHATTERSEEK_CHROMS, , drop = FALSE]
      raw <- raw[order(match(raw$chrom, SHATTERSEEK_CHROMS), as.numeric(raw$start)), , drop = FALSE]
      pick <- function(nm, default = NA) if (nm %in% colnames(raw)) suppressWarnings(as.numeric(raw[[nm]])) else rep(default, nrow(raw))
      if (all(c("nMaj1_A", "nMin1_A") %in% colnames(raw))) {
        nMaj1 <- pick("nMaj1_A"); nMin1 <- pick("nMin1_A"); frac1 <- pick("frac1_A", 1)
        nMaj2 <- pick("nMaj2_A"); nMin2 <- pick("nMin2_A"); frac2 <- pick("frac2_A")
        SD <- pick("SDfrac_A")
      } else if (all(c("nMaj1", "nMin1") %in% colnames(raw))) {
        nMaj1 <- pick("nMaj1"); nMin1 <- pick("nMin1"); frac1 <- pick("frac1", 1)
        nMaj2 <- pick("nMaj2"); nMin2 <- pick("nMin2"); frac2 <- pick("frac2"); SD <- pick("SD")
      } else {
        cn <- read_cn_file(m[[cn_col]][i], filter_pass = filter_pass)   # total CN only
        cn$chrom <- normalize_chr_bare(cn$chrom); cn <- cn[cn$chrom %in% SHATTERSEEK_CHROMS, ]
        cn <- cn[order(match(cn$chrom, SHATTERSEEK_CHROMS), cn$start), ]
        raw <- data.frame(chrom = cn$chrom, start = cn$start, end = cn$end)
        tot <- round(cn$total_cn); nMaj1 <- ceiling(tot / 2); nMin1 <- floor(tot / 2); frac1 <- rep(1, nrow(raw))
        nMaj2 <- nMin2 <- frac2 <- SD <- rep(NA, nrow(raw))
      }
      frac1[is.na(frac1)] <- 1
      tot <- nMaj1 + nMin1
      tot2 <- nMaj2 + nMin2
      tot_w <- ifelse(is.na(tot2) | is.na(frac2), tot, frac1 * tot + (1 - frac1) * tot2)
      ploidy <- if ("ploidy" %in% colnames(raw)) as.numeric(raw$ploidy[1]) else {
        auto <- raw$chrom %in% as.character(1:22) & !is.na(tot_w)
        len <- as.numeric(raw$end[auto]) - as.numeric(raw$start[auto]) + 1
        round(sum(tot_w[auto] * len) / sum(len), 3)
      }
      cn_rows[[sid]] <- data.frame(samplename = sid, chr = raw$chrom, start = as.numeric(raw$start), end = as.numeric(raw$end),
                                   nMaj1 = nMaj1, nMin1 = nMin1, frac1 = frac1, nMaj2 = nMaj2, nMin2 = nMin2, frac2 = frac2,
                                   SD = SD, ploidy = ploidy, stringsAsFactors = FALSE)
      sprintf("%d SVs, %d CN segments, ploidy %.2f", nrow(sv), nrow(raw), ploidy)
    }, silent = TRUE)
    status <- if (inherits(res, "try-error")) paste("FAILED:", conditionMessage(attr(res, "condition"))) else res
    log[[sid]] <- data.frame(sampleID = sid, status = status, stringsAsFactors = FALSE)
    message(sprintf("[%d/%d] %s: %s", i, nrow(m), sid, status))
  }
  sv_all <- do.call(rbind, sv_rows); cn_all <- do.call(rbind, cn_rows)
  write.table(sv_all, file.path(outdir, "SVs.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(cn_all, file.path(outdir, "CNAs.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(do.call(rbind, log), file.path(outdir, "build_cohort_tables_log.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote ", file.path(outdir, "SVs.txt"), " (", nrow(sv_all), " SVs) and ", file.path(outdir, "CNAs.txt"),
          " (", nrow(cn_all), " segments) for ", length(cn_rows), " sample(s)")
}

main()
