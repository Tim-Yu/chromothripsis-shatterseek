# ---------------------------------------------------------------------------
# shatterseek_lib.R
#
# Shared helpers for running ShatterSeek (Cortes-Ciriano et al. 2020,
# https://github.com/parklab/ShatterSeek) on the group's SV / copy-number data.
#
# Sourced by
#   run_shatterseek_persample.R  - per-sample bedpe + Battenberg/DRAGEN CN files
#                                  (the "circos pipeline" layout)
#   run_shatterseek_cohort.R     - cohort tables SVs.txt / CNAs.txt
#                                  (the legacy chromothripsis.r layout)
#
# Everything here is plain R; no side effects at source time.
# ---------------------------------------------------------------------------

# Library paths for the GEL research environment, applied only when the
# directories exist, so nothing changes elsewhere. Equivalent to
#   .libPaths(c("/home/byu/R/x86_64-pc-linux-gnu-library/4.5", .libPaths(),
#               "/tools/aws-workspace-ubuntu-apps/ce/R/4.5.3", "/tools/aws-workspace-apps/ce/R/4.2.1/"))
# i.e. the personal library holding ShatterSeek first, the shared trees with the
# dependencies last. Override / extend with SHATTERSEEK_LIB (prepended) and
# SHATTERSEEK_EXTRA_LIBS="dir1:dir2" (appended).
local({
  personal <- c(Sys.getenv("SHATTERSEEK_LIB", ""), "/home/byu/R/x86_64-pc-linux-gnu-library/4.5")
  extra <- c("/tools/aws-workspace-ubuntu-apps/ce/R/4.5.3", "/tools/aws-workspace-apps/ce/R/4.2.1/",
             strsplit(Sys.getenv("SHATTERSEEK_EXTRA_LIBS", ""), ":")[[1]])
  personal <- personal[nzchar(personal) & dir.exists(personal)]
  extra <- extra[nzchar(extra) & dir.exists(extra)]
  if (length(personal) || length(extra)) .libPaths(c(personal, .libPaths(), extra))
})

suppressPackageStartupMessages({
  library(ShatterSeek)
  library(GenomicRanges)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x[1]) || x[1] == "") y else x

parse_bool <- function(x, default = FALSE) {
  if (is.null(x) || is.na(x) || x == "") return(default)
  tolower(as.character(x)) %in% c("true", "1", "yes", "y")
}

# Generic "--key value" parser shared by all scripts (flags listed in `flags`
# take no value).
parse_cli <- function(flags = c("help")) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Invalid argument: ", key)
    key <- sub("^--", "", key)
    if (key %in% flags) { out[[key]] <- TRUE; i <- i + 1; next }
    if (i == length(args)) stop("Missing value for argument --", key)
    out[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

# ---------------------------------------------------------------------------
# Genome constants
# ---------------------------------------------------------------------------
SHATTERSEEK_CHROMS <- c(as.character(1:22), "X")

# ShatterSeek wants bare chromosome names ("1".."22","X"); 23 -> X, 24 -> Y.
normalize_chr_bare <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x[x == "23"] <- "X"
  x[x == "24"] <- "Y"
  x
}

# ---------------------------------------------------------------------------
# SV class / strand harmonisation
#
# ShatterSeek SVtype must be one of DEL, DUP, h2hINV, t2tINV, TRA and the
# strands must be "+"/"-".  Our consensus bedpe files carry svclass in that
# vocabulary already (sometimes "INV" or "BND"), and strands wrapped in
# brackets ("[+]").  The legacy cohort tables carry bare +/- strands.
#
# Strand convention used by ShatterSeek / Korbel & Campbell for intrachromosomal
# joins:  DEL = +/-,  DUP = -/+,  h2hINV = +/+,  t2tINV = -/-.
# When svclass is informative we derive strands from it (this is also what
# makes the fragment-joins test meaningful); when it is not (TRA, INV, BND)
# we keep the strands given in the file.
# ---------------------------------------------------------------------------
CLASS_STRANDS <- list(DEL = c("+", "-"), DUP = c("-", "+"), h2hINV = c("+", "+"), t2tINV = c("-", "-"))

clean_strand <- function(x) {
  x <- gsub("[^+-]", "", as.character(x))   # strips "[" "]" etc.
  x[!x %in% c("+", "-")] <- NA
  x
}

canonical_svtype <- function(svclass, strand1, strand2, same_chr) {
  s <- tolower(trimws(as.character(svclass)))
  out <- rep(NA_character_, length(s))
  out[s %in% c("del", "deletion")] <- "DEL"
  out[s %in% c("dup", "duplication", "tandem_dup", "tdup")] <- "DUP"
  out[s %in% c("h2hinv", "h2h", "inv_h2h")] <- "h2hINV"
  out[s %in% c("t2tinv", "t2t", "inv_t2t")] <- "t2tINV"
  out[s %in% c("tra", "bnd", "ctx", "translocation", "itx")] <- "TRA"
  # ambiguous "INV": resolve by strand if possible, else assign h2hINV
  inv <- s %in% c("inv", "inversion")
  out[inv & strand1 == "+" & strand2 == "+"] <- "h2hINV"
  out[inv & strand1 == "-" & strand2 == "-"] <- "t2tINV"
  out[inv & is.na(out)] <- "h2hINV"
  # unknown class but strands available: infer from strands
  unk <- is.na(out) & !is.na(strand1) & !is.na(strand2)
  key <- paste0(strand1, strand2)
  out[unk & key == "+-"] <- "DEL"
  out[unk & key == "-+"] <- "DUP"
  out[unk & key == "++"] <- "h2hINV"
  out[unk & key == "--"] <- "t2tINV"
  # anything interchromosomal is a TRA for ShatterSeek regardless of label
  out[!same_chr] <- "TRA"
  out
}

# Build the ShatterSeek SVs object from a harmonised data.frame with columns
# chrom1,pos1,chrom2,pos2,svclass,strand1,strand2 (strands may be NA/garbage).
build_sv_object <- function(df, keep_chroms = SHATTERSEEK_CHROMS, derive_strands_from_class = TRUE) {
  df$chrom1 <- normalize_chr_bare(df$chrom1)
  df$chrom2 <- normalize_chr_bare(df$chrom2)
  df$pos1 <- as.numeric(df$pos1)
  df$pos2 <- as.numeric(df$pos2)
  df$strand1 <- clean_strand(df$strand1)
  df$strand2 <- clean_strand(df$strand2)
  same <- df$chrom1 == df$chrom2
  df$SVtype <- canonical_svtype(df$svclass, df$strand1, df$strand2, same)

  if (derive_strands_from_class) {
    for (k in names(CLASS_STRANDS)) {
      idx <- which(df$SVtype == k)
      df$strand1[idx] <- CLASS_STRANDS[[k]][1]
      df$strand2[idx] <- CLASS_STRANDS[[k]][2]
    }
  }
  # TRA with missing strands: fill at random so the object validates
  # (strand of TRA only feeds the inter-chromosomal fragment-join count).
  miss <- is.na(df$strand1); if (any(miss)) df$strand1[miss] <- sample(c("+", "-"), sum(miss), TRUE)
  miss <- is.na(df$strand2); if (any(miss)) df$strand2[miss] <- sample(c("+", "-"), sum(miss), TRUE)

  # order intrachromosomal breakpoints so pos1 <= pos2 (cluster.SV builds IRanges(pos1,pos2))
  swap <- same & df$pos1 > df$pos2
  if (any(swap)) {
    tmp <- df$pos1[swap]; df$pos1[swap] <- df$pos2[swap]; df$pos2[swap] <- tmp
    tmp <- df$strand1[swap]; df$strand1[swap] <- df$strand2[swap]; df$strand2[swap] <- tmp
  }

  keep <- df$chrom1 %in% keep_chroms & df$chrom2 %in% keep_chroms &
    !is.na(df$pos1) & !is.na(df$pos2) & !is.na(df$SVtype)
  dropped <- sum(!keep)
  df <- df[keep, , drop = FALSE]
  # ShatterSeek chokes on zero-length intrachromosomal ranges
  zero <- df$chrom1 == df$chrom2 & df$pos1 == df$pos2
  df <- df[!zero, , drop = FALSE]
  df <- unique(df[, c("chrom1", "pos1", "chrom2", "pos2", "SVtype", "strand1", "strand2")])

  obj <- SVs(chrom1 = as.character(df$chrom1), pos1 = as.numeric(df$pos1),
             chrom2 = as.character(df$chrom2), pos2 = as.numeric(df$pos2),
             SVtype = as.character(df$SVtype),
             strand1 = as.character(df$strand1), strand2 = as.character(df$strand2))
  attr(obj, "n_dropped") <- dropped
  attr(obj, "df") <- df
  obj
}

# Build the CNVsegs object from a harmonised data.frame chrom,start,end,total_cn
build_cn_object <- function(df, keep_chroms = SHATTERSEEK_CHROMS, round_cn = TRUE) {
  df$chrom <- normalize_chr_bare(df$chrom)
  df$start <- as.numeric(df$start)
  df$end <- as.numeric(df$end)
  df$total_cn <- as.numeric(df$total_cn)
  df <- df[df$chrom %in% keep_chroms & !is.na(df$total_cn) & !is.na(df$start) & !is.na(df$end) & df$end >= df$start, ]
  if (round_cn) df$total_cn <- round(df$total_cn)
  df <- df[order(match(df$chrom, keep_chroms), df$start), ]
  # merge adjacent segments with identical CN (segmentation artefacts inflate
  # "number of segments" and break the oscillation count)
  df <- merge_equal_adjacent(df)
  CNVsegs(chrom = as.character(df$chrom), start = df$start, end = df$end, total_cn = df$total_cn)
}

merge_equal_adjacent <- function(df) {
  if (nrow(df) < 2) return(df)
  keep <- rep(TRUE, nrow(df))
  for (i in 2:nrow(df)) {
    j <- max(which(keep[1:(i - 1)]))
    if (df$chrom[i] == df$chrom[j] && df$total_cn[i] == df$total_cn[j]) {
      df$end[j] <- df$end[i]
      keep[i] <- FALSE
    }
  }
  df[keep, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Readers for the per-sample "circos pipeline" formats
# ---------------------------------------------------------------------------

# Consensus-caller bedpe: chrom1 start1 end1 chrom2 start2 end2 sv_id pe_support
# strand1 strand2 svclass svmethod  (header may be absent -> positional)
read_sv_bedpe <- function(path) {
  first <- readLines(path, n = 1)
  has_header <- grepl("chrom1|CHROM_A|start1", first, ignore.case = TRUE)
  sv <- read.delim(path, header = has_header, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
  if (!has_header) {
    nm <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "sv_id", "pe_support", "strand1", "strand2", "svclass", "svmethod")
    colnames(sv)[seq_len(min(ncol(sv), length(nm)))] <- nm[seq_len(min(ncol(sv), length(nm)))]
  }
  # VCF2BEDPE-style capitalised headers
  alias <- c(CHROM_A = "chrom1", START_A = "start1", END_A = "end1", CHROM_B = "chrom2", START_B = "start2", END_B = "end2",
             STRAND_A = "strand1", STRAND_B = "strand2", SVTYPE = "svtype_raw", ID = "sv_id")
  for (a in names(alias)) if (a %in% colnames(sv) && !alias[[a]] %in% colnames(sv)) colnames(sv)[colnames(sv) == a] <- alias[[a]]
  if (!"svclass" %in% colnames(sv)) sv$svclass <- if ("svtype_raw" %in% colnames(sv)) sv$svtype_raw else NA
  req <- c("chrom1", "start1", "chrom2", "start2")
  miss <- setdiff(req, colnames(sv))
  if (length(miss)) stop("SV file ", path, " missing columns: ", paste(miss, collapse = ", "))
  if (!"strand1" %in% colnames(sv)) sv$strand1 <- NA
  if (!"strand2" %in% colnames(sv)) sv$strand2 <- NA
  data.frame(chrom1 = sv$chrom1, pos1 = as.numeric(sv$start1), chrom2 = sv$chrom2, pos2 = as.numeric(sv$start2),
             svclass = sv$svclass, strand1 = sv$strand1, strand2 = sv$strand2, stringsAsFactors = FALSE)
}

# Copy number: Battenberg subclones / gapfilled file (chr startpos endpos ... ntot|tot_n, nMaj1_A nMin1_A frac1_A nMaj2_A ...)
# or DRAGEN-style export (Chromosome, Start Position, End Position, Tumour TCN, vcf_filter).
CN_COL_ALIASES <- list(
  chrom = c("chrom", "chr", "Chromosome", "chromosome", "seqnames"),
  start = c("start", "startpos", "Start Position", "Start", "start_pos"),
  end = c("end", "endpos", "End Position", "End", "end_pos"),
  total_cn = c("total_cn", "ntot", "tot_n", "Tumour TCN", "TCN", "nTot", "total_copy_number", "cn")
)

resolve_cn_cols <- function(cn) {
  for (canon in names(CN_COL_ALIASES)) {
    if (canon %in% colnames(cn)) next
    hit <- intersect(CN_COL_ALIASES[[canon]], colnames(cn))
    if (length(hit)) colnames(cn)[colnames(cn) == hit[[1]]] <- canon
  }
  cn
}

read_cn_file <- function(path, filter_pass = FALSE, subclonal = c("clonal", "weighted")) {
  subclonal <- match.arg(subclonal)
  cn <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  cn <- resolve_cn_cols(cn)
  if (filter_pass && "vcf_filter" %in% colnames(cn)) cn <- cn[tolower(cn$vcf_filter) == "pass", , drop = FALSE]
  if (!"total_cn" %in% colnames(cn)) {
    # Battenberg subclones without an ntot column: use the clonal (or
    # fraction-weighted) total of major + minor allele copy numbers.
    if (all(c("nMaj1_A", "nMin1_A") %in% colnames(cn))) {
      tot1 <- cn$nMaj1_A + cn$nMin1_A
      if (subclonal == "weighted" && all(c("nMaj2_A", "nMin2_A", "frac1_A") %in% colnames(cn))) {
        tot2 <- cn$nMaj2_A + cn$nMin2_A
        f1 <- ifelse(is.na(cn$frac1_A), 1, cn$frac1_A)
        cn$total_cn <- ifelse(is.na(tot2), tot1, f1 * tot1 + (1 - f1) * tot2)
      } else cn$total_cn <- tot1
    } else if (all(c("nMaj1", "nMin1") %in% colnames(cn))) {
      cn$total_cn <- cn$nMaj1 + cn$nMin1
    } else if ("Seg.CN" %in% colnames(cn)) {
      cn$total_cn <- 2 * 2^cn$`Seg.CN`      # log2 ratio -> absolute
    } else stop("CN file ", path, " has no usable total copy-number column")
  }
  miss <- setdiff(c("chrom", "start", "end", "total_cn"), colnames(cn))
  if (length(miss)) stop("CN file ", path, " missing columns: ", paste(miss, collapse = ", "))
  data.frame(chrom = cn$chrom, start = as.numeric(cn$start), end = as.numeric(cn$end), total_cn = as.numeric(cn$total_cn),
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Readers for the legacy cohort tables
#   SVs.txt : samplename chr1 strand1 pos1 chr2 strand2 pos2 svclass
#   CNAs.txt: samplename chr start end nMaj1 nMin1 frac1 nMaj2 nMin2 frac2 SD ploidy
# (column *positions* are what the legacy script relied on; we accept either
#  the names above or positional order)
# ---------------------------------------------------------------------------
read_cohort_sv <- function(path) {
  sv <- read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  nm <- c("samplename", "chr1", "strand1", "pos1", "chr2", "strand2", "pos2", "svclass")
  if (!all(nm %in% colnames(sv))) colnames(sv)[seq_along(nm)] <- nm
  sv
}

read_cohort_cn <- function(path, subclonal = c("clonal", "weighted")) {
  subclonal <- match.arg(subclonal)
  cn <- read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  nm <- c("samplename", "chr", "start", "end", "nMaj1", "nMin1", "frac1", "nMaj2", "nMin2", "frac2", "SD", "ploidy")
  if (!all(nm[1:6] %in% colnames(cn))) colnames(cn)[seq_len(min(ncol(cn), length(nm)))] <- nm[seq_len(min(ncol(cn), length(nm)))]
  tot1 <- as.numeric(cn$nMaj1) + as.numeric(cn$nMin1)
  if (subclonal == "weighted" && all(c("nMaj2", "nMin2", "frac1") %in% colnames(cn))) {
    tot2 <- as.numeric(cn$nMaj2) + as.numeric(cn$nMin2)
    f1 <- ifelse(is.na(cn$frac1), 1, as.numeric(cn$frac1))
    cn$total_cn <- ifelse(is.na(tot2), tot1, f1 * tot1 + (1 - f1) * tot2)
  } else cn$total_cn <- tot1
  cn
}

# ---------------------------------------------------------------------------
# Runtime patch for ShatterSeek 1.1 statistical_criteria():
#   when collecting inter-chromosomal SVs whose *second* breakpoint lies on the
#   candidate chromosome, the code tests inter$pos1 (the other chromosome's
#   position) instead of inter$pos2. Translocations recorded as
#   chrom1=partner, chrom2=candidate are therefore mostly dropped, which
#   under-counts number_TRA / clusterSize_including_TRA and weakens the HC2
#   criterion for the candidate chromosome. We fix the comparison in place.
# Verified against the function text; the patch is skipped (with a warning)
# if the expected line is not found exactly once.
# ---------------------------------------------------------------------------
patch_shatterseek_inter_bug <- function(verbose = TRUE) {
  if (isTRUE(getOption("shatterseek.patched"))) return(invisible(TRUE))
  f <- get("statistical_criteria", envir = asNamespace("ShatterSeek"))
  txt <- deparse(f, width.cutoff = 500)
  pat <- "idx_inter2 = which(inter$chrom2 == cand & inter$pos1 >= (min_now - 10000) & inter$pos2 <= (max_now + 10000))"
  hits <- grep(pat, txt, fixed = TRUE)
  if (length(hits) != 1) { warning("ShatterSeek inter-chromosomal patch not applied (pattern found ", length(hits), " times)"); return(invisible(FALSE)) }
  txt[hits] <- sub(pat, "idx_inter2 = which(inter$chrom2 == cand & inter$pos2 >= (min_now - 10000) & inter$pos2 <= (max_now + 10000))", txt[hits], fixed = TRUE)
  g <- eval(parse(text = txt), envir = asNamespace("ShatterSeek"))
  environment(g) <- asNamespace("ShatterSeek")
  utils::assignInNamespace("statistical_criteria", g, ns = "ShatterSeek")
  options(shatterseek.patched = TRUE)
  if (verbose) message("ShatterSeek: applied inter-chromosomal pos1/pos2 patch to statistical_criteria()")
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Baseline-return oscillation (DFSP ring / amplicon pattern)
#
# ShatterSeek counts oscillations between exactly 2 (or 3 adjacent) CN states.
# DFSP ring chromosomes and other high-level amplicons alternate between ONE
# baseline state and VARYING amplified states (e.g. 3,11,3,9,3,11,3,10,3), so
# the 2-state count stays small although the profile is clearly oscillating.
# Here we report the longest run of consecutive segments in the cluster region
# in which every other segment returns to the same CN while the segments in
# between differ from it (by >= 1 copy). The 2-state oscillation is a special
# case of this, so baseline_osc >= 2-state osc always.
# ---------------------------------------------------------------------------
baseline_oscillation <- function(cn) {
  n <- length(cn)
  if (n < 3) return(n)
  best <- 0
  for (s in 1:(n - 2)) {
    base <- cn[s]
    len <- 1
    i <- s
    while (i + 2 <= n && cn[i + 2] == base && cn[i + 1] != base) { len <- len + 2; i <- i + 2 }
    if (len > best) best <- len
  }
  max(best, 1)
}

# CN segments overlapping a cluster region on one chromosome (same selection as ShatterSeek)
region_cn <- function(cn_df, chrom, start, end) {
  x <- cn_df[cn_df$chrom == chrom, ]
  x <- x[order(x$start), ]
  x[x$end >= start & x$start <= end, , drop = FALSE]
}

add_baseline_oscillation <- function(summary_df, cn_df) {
  summary_df$max_number_oscillating_CN_segments_baseline <- NA_real_
  summary_df$region_CN_states <- NA_character_
  for (i in seq_len(nrow(summary_df))) {
    if (is.na(summary_df$start[i])) next
    r <- region_cn(cn_df, summary_df$chrom[i], summary_df$start[i], summary_df$end[i])
    if (nrow(r) == 0) next
    summary_df$max_number_oscillating_CN_segments_baseline[i] <- baseline_oscillation(r$total_cn)
    summary_df$region_CN_states[i] <- paste(r$total_cn, collapse = ";")
  }
  summary_df
}

# ---------------------------------------------------------------------------
# Running ShatterSeek on one sample
# ---------------------------------------------------------------------------
run_shatterseek_sample <- function(sv_df, cn_df, genome = "hg38", min_size = 1, round_cn = TRUE, seed = 1, patch = TRUE) {
  if (patch) patch_shatterseek_inter_bug(verbose = FALSE)
  sv_obj <- build_sv_object(sv_df)
  cn_obj <- build_cn_object(cn_df, round_cn = round_cn)
  if (length(sv_obj@chrom1) == 0) stop("no usable SVs after filtering")
  # the 'random distribution of breakpoints' test in ShatterSeek draws a random
  # exponential sample, so fix the seed for reproducible p-values
  set.seed(seed)
  res <- suppressMessages(suppressWarnings(
    utils::capture.output(out <- shatterseek(SV.sample = sv_obj, seg.sample = cn_obj, min.Size = min_size, genome = genome))
  ))
  attr(out, "sv_df") <- attr(sv_obj, "df")
  attr(out, "cn_df") <- as(cn_obj, "data.frame")
  attr(out, "n_sv_dropped") <- attr(sv_obj, "n_dropped")
  out
}

# ---------------------------------------------------------------------------
# Calling rules
#
# Published ShatterSeek criteria (Cortes-Ciriano et al. 2020, Nat Genet;
# ShatterSeek tutorial section "Criteria to call chromothripsis"):
#   HC1: >=6 interleaved intra-chr SVs, >=7 oscillating CN segments (2 states),
#        fragment-joins test NOT rejected (p > 0.05), and (chr breakpoint
#        enrichment p < 0.05 OR exponential breakpoint distribution p < 0.05)
#   HC2: >=3 interleaved intra-chr SVs AND >=4 inter-chr SVs, >=7 oscillating
#        CN segments (2 states), fragment-joins test NOT rejected
#   LC : >=6 interleaved intra-chr SVs, 4-6 oscillating CN segments (2 states),
#        fragment-joins NOT rejected, and (enrichment OR exponential test)
#
# Extension for DFSP-type ring / high-level amplicon events (the GMS cases):
# these oscillate across MANY CN states (3, 5, 9, 11...), so the strict 2-state
# count under-reports. We therefore also report an "Extended" tier using the
# 3-state oscillation count and the SV cluster size including translocations,
# which is flagged separately so it can be reviewed and not confused with the
# published calls. All thresholds are parameters so they can be tuned.
# ---------------------------------------------------------------------------
default_thresholds <- function() list(
  hc_min_intra = 6, hc_min_intra_with_tra = 3, hc_min_tra = 4, hc_min_osc2 = 7,
  lc_min_intra = 6, lc_min_osc2 = 4, lc_max_osc2 = 6,
  p_joins = 0.05, p_enrich = 0.05, p_exp = 0.05,
  ext_min_cluster_incl_tra = 10, ext_min_osc3 = 7, ext_min_osc_baseline = 7, ext_min_cn_segments = 10, ext_min_tra = 4,
  use_3state_for_lc = FALSE
)

classify_summary <- function(s, th = default_thresholds()) {
  g <- function(x) ifelse(is.na(x), NA, x)
  intra <- s$number_DEL + s$number_DUP + s$number_h2hINV + s$number_t2tINV
  tra <- g(s$number_TRA); tra[is.na(tra)] <- 0
  osc2 <- g(s$max_number_oscillating_CN_segments_2_states)
  osc3 <- g(s$max_number_oscillating_CN_segments_3_states)
  joins_ok <- is.na(s$pval_fragment_joins) | s$pval_fragment_joins > th$p_joins
  # when fragment-joins p is NA (e.g. no classified intra SVs) treat as not-rejected
  enrich_ok <- !is.na(s$chr_breakpoint_enrichment) & s$chr_breakpoint_enrichment < th$p_enrich
  exp_ok <- (!is.na(s$pval_exp_cluster) & s$pval_exp_cluster < th$p_exp) | (!is.na(s$pval_exp_chr) & s$pval_exp_chr < th$p_exp)
  osc2_ok_hc <- !is.na(osc2) & osc2 >= th$hc_min_osc2
  osc_lc <- if (isTRUE(th$use_3state_for_lc)) osc3 else osc2
  osc_ok_lc <- !is.na(osc_lc) & osc_lc >= th$lc_min_osc2 & osc_lc <= th$lc_max_osc2

  hc1 <- intra >= th$hc_min_intra & osc2_ok_hc & joins_ok & (enrich_ok | exp_ok)
  hc2 <- intra >= th$hc_min_intra_with_tra & tra >= th$hc_min_tra & osc2_ok_hc & joins_ok
  lc <- intra >= th$lc_min_intra & osc_ok_lc & joins_ok & (enrich_ok | exp_ok)
  oscb <- if ("max_number_oscillating_CN_segments_baseline" %in% colnames(s)) g(s$max_number_oscillating_CN_segments_baseline) else rep(NA, nrow(s))
  ext <- !is.na(s$clusterSize_including_TRA) & s$clusterSize_including_TRA >= th$ext_min_cluster_incl_tra &
    ((!is.na(osc3) & osc3 >= th$ext_min_osc3) | (!is.na(oscb) & oscb >= th$ext_min_osc_baseline) |
       (!is.na(s$number_CNV_segments) & s$number_CNV_segments >= th$ext_min_cn_segments & tra >= th$ext_min_tra)) &
    (enrich_ok | exp_ok | tra >= th$ext_min_tra)

  call <- rep("None", nrow(s))
  call[ext] <- "Extended_ringlike"
  call[lc] <- "LowConfidence"
  call[hc1 | hc2] <- "HighConfidence"
  hc1[is.na(hc1)] <- FALSE; hc2[is.na(hc2)] <- FALSE; lc[is.na(lc)] <- FALSE; ext[is.na(ext)] <- FALSE
  data.frame(number_intra_SVs = intra, criterion_HC1 = hc1, criterion_HC2 = hc2, criterion_LC = lc,
             criterion_extended = ext, call = call, stringsAsFactors = FALSE)
}

# Tidy per-chromosome table for one sample
summarise_sample <- function(out, sample_id, patient_id = NA, th = default_thresholds()) {
  s <- out@chromSummary
  if (!is.null(attr(out, "cn_df"))) s <- add_baseline_oscillation(s, attr(out, "cn_df"))
  cls <- classify_summary(s, th)
  res <- cbind(data.frame(patientID = patient_id, sampleID = sample_id, stringsAsFactors = FALSE), s, cls)
  res$region <- ifelse(is.na(res$start), NA, paste0("chr", res$chrom, ":", res$start, "-", res$end))
  res
}

# ---------------------------------------------------------------------------
# Plotting: ShatterSeek's own chromosome plot (SV arcs + CN) for called chromosomes
# ---------------------------------------------------------------------------
plot_called_chromosomes <- function(out, summary_df, outdir, sample_id, genome = "hg38",
                                    which_calls = c("HighConfidence", "LowConfidence", "Extended_ringlike"),
                                    width = 10, height = 6) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  chrs <- summary_df$chrom[summary_df$call %in% which_calls]
  files <- character(0)
  for (chr in chrs) {
    f <- file.path(outdir, paste0(sample_id, "_chr", chr, "_shatterseek.pdf"))
    ok <- try({
      p <- plot_chromothripsis(ShatterSeek_output = out, chr = as.character(chr), sample_name = sample_id, genome = genome)
      pdf(f, width = width, height = height)
      grid::grid.draw(gridExtra::arrangeGrob(p[[1]], p[[2]], p[[3]], p[[4]], nrow = 4, ncol = 1, heights = c(0.2, 0.4, 0.4, 0.4)))
      dev.off()
    }, silent = TRUE)
    if (inherits(ok, "try-error")) { try(dev.off(), silent = TRUE); warning("plot failed for ", sample_id, " chr", chr, ": ", ok) } else files <- c(files, f)
  }
  invisible(files)
}

# Write the standard output set for one sample
write_sample_outputs <- function(out, summary_df, outdir, sample_id, save_rds = TRUE) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  write.table(summary_df, file.path(outdir, paste0(sample_id, "_shatterseek_summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  if (save_rds) saveRDS(out, file.path(outdir, paste0(sample_id, "_shatterseek.rds")))
  invisible(NULL)
}

# Columns kept in the cohort-level "calls" table (one row per chromosome with a call)
CALL_COLUMNS <- c("patientID", "sampleID", "chrom", "region", "start", "end", "call",
                  "criterion_HC1", "criterion_HC2", "criterion_LC", "criterion_extended",
                  "number_intra_SVs", "number_DEL", "number_DUP", "number_h2hINV", "number_t2tINV", "number_TRA",
                  "clusterSize_including_TRA", "number_SVs_sample", "number_CNV_segments",
                  "max_number_oscillating_CN_segments_2_states", "max_number_oscillating_CN_segments_3_states",
                  "max_number_oscillating_CN_segments_baseline", "region_CN_states",
                  "pval_fragment_joins", "chr_breakpoint_enrichment", "pval_exp_chr", "pval_exp_cluster",
                  "inter_other_chroms", "inter_other_chroms_coords_all")

# ---------------------------------------------------------------------------
# Threshold overrides from the command line: "--thresholds hc_min_osc2=6,p_joins=0.1"
# ---------------------------------------------------------------------------
parse_thresholds <- function(spec, th = default_thresholds()) {
  if (is.null(spec) || spec == "") return(th)
  for (kv in strsplit(spec, ",")[[1]]) {
    p <- strsplit(trimws(kv), "=")[[1]]
    if (length(p) != 2) stop("Bad --thresholds entry: ", kv)
    if (!p[1] %in% names(th)) stop("Unknown threshold: ", p[1], ". Known: ", paste(names(th), collapse = ", "))
    th[[p[1]]] <- if (is.logical(th[[p[1]]])) parse_bool(p[2]) else as.numeric(p[2])
  }
  th
}

# ---------------------------------------------------------------------------
# Driver shared by both runners.
#   samples: data.frame with patientID, sampleID and a loader per row supplied
#            through `load_fun(i)` returning list(sv = sv_df, cn = cn_df)
# Writes per-sample outputs and the two cohort tables; returns the full table.
# ---------------------------------------------------------------------------
run_shatterseek_batch <- function(samples, load_fun, outdir, genome = "hg38", th = default_thresholds(),
                                  min_size = 1, do_plot = TRUE, save_rds = TRUE, plot_calls = c("HighConfidence", "LowConfidence", "Extended_ringlike")) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  per_dir <- file.path(outdir, "per_sample")
  plot_dir <- file.path(outdir, "plots")
  all_rows <- list(); log_rows <- list()
  for (i in seq_len(nrow(samples))) {
    sid <- samples$sampleID[i]; pid <- samples$patientID[i]
    message(sprintf("[%d/%d] %s / %s", i, nrow(samples), pid, sid))
    status <- "ok"; msg <- ""
    res <- try({
      d <- load_fun(i)
      out <- run_shatterseek_sample(d$sv, d$cn, genome = genome, min_size = min_size)
      summ <- summarise_sample(out, sid, pid, th)
      write_sample_outputs(out, summ, per_dir, sid, save_rds = save_rds)
      if (do_plot) plot_called_chromosomes(out, summ, plot_dir, sid, genome = genome, which_calls = plot_calls)
      all_rows[[sid]] <- summ
      msg <- sprintf("%d SVs used, %d dropped; calls: %s", length(out@detail$SV$chrom1) + nrow(out@detail$SVinter),
                     attr(out, "n_sv_dropped"), paste(summ$call[summ$call != "None"], collapse = ","))
    }, silent = TRUE)
    if (inherits(res, "try-error")) { status <- "failed"; msg <- conditionMessage(attr(res, "condition")); warning(sid, ": ", msg) }
    log_rows[[sid]] <- data.frame(patientID = pid, sampleID = sid, status = status, message = msg, stringsAsFactors = FALSE)
    message("    ", status, ": ", msg)
  }
  all_df <- if (length(all_rows)) do.call(rbind, all_rows) else NULL
  rownames(all_df) <- NULL
  write.table(do.call(rbind, log_rows), file.path(outdir, "shatterseek_run_log.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  if (!is.null(all_df)) {
    write.table(all_df, file.path(outdir, "shatterseek_all_chromosomes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    calls <- all_df[all_df$call != "None", intersect(CALL_COLUMNS, colnames(all_df)), drop = FALSE]
    write.table(calls, file.path(outdir, "shatterseek_calls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    message("\nWrote: ", file.path(outdir, "shatterseek_calls.tsv"), " (", nrow(calls), " called chromosome(s) in ", length(all_rows), " sample(s))")
  }
  write_thresholds(th, file.path(outdir, "thresholds_used.txt"))
  invisible(all_df)
}

write_thresholds <- function(th, path) {
  writeLines(paste0(names(th), "=", unlist(th)), path)
}
