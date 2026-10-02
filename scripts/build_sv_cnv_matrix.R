#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || is.na(x) || x == "") y else x
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()

  if (length(args) == 0) {
    return(out)
  }

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Invalid argument: ", key)
    }
    key <- sub("^--", "", key)

    if (key == "help") {
      out$help <- TRUE
      i <- i + 1
      next
    }

    if (i == length(args)) {
      stop("Missing value for argument --", key)
    }

    out[[key]] <- args[[i + 1]]
    i <- i + 2
  }

  out
}

print_usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript build_sv_cnv_matrix.R \\\n",
    "    --patients <patients_txt> \\\n",
    "    --sv-root <sv_root_dir> \\\n",
    "    --cn-root <cn_root_dir> \\\n",
    "    --out <output_tsv> \\\n",
    "    [--sv-file-pattern <regex>] [--cn-dir-pattern <regex>] [--cn-file-pattern <regex>]\n\n",
    "Notes:\n",
    "  - patients file: one patient ID per line\n",
    "  - SV search: <sv_root>/<patientID>/<sampleID>/<sv-file-pattern match>\n",
    "      default --sv-file-pattern: {sampleID}\\.somatic\\.sv\\.bedpe\n",
    "  - CN dir search (first level under --cn-root): <cn-dir-pattern match>\n",
    "      default --cn-dir-pattern: P_{patientID}_T_{sampleID}_.*\n",
    "  - CN file search inside matched CN dir: <cn-file-pattern match>\n",
    "      default --cn-file-pattern: {sampleID}_copynumber_1_gapfilled_LogR_new\\.txt\n",
    "  - Patterns are regex (anchored with ^ and $ automatically) and support\n",
    "    placeholders {patientID} and {sampleID}\n",
    sep = ""
  )
}

fill_placeholders <- function(pattern, patient_id, sample_id) {
  pattern <- gsub("{patientID}", patient_id, pattern, fixed = TRUE)
  pattern <- gsub("{sampleID}", sample_id, pattern, fixed = TRUE)
  pattern
}

find_matching_entries <- function(base_dir, pattern, want_dir) {
  entries <- list.files(base_dir, full.names = TRUE, recursive = FALSE, all.files = FALSE, no.. = TRUE)
  entries <- entries[file.info(entries)$isdir %in% want_dir]

  anchored <- pattern
  if (!startsWith(anchored, "^")) anchored <- paste0("^", anchored)
  if (!endsWith(anchored, "$")) anchored <- paste0(anchored, "$")

  sort(entries[grepl(anchored, basename(entries))])
}

read_patients <- function(path) {
  if (!file.exists(path)) {
    stop("Patients file not found: ", path)
  }
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}

list_sample_dirs <- function(patient_dir) {
  if (!dir.exists(patient_dir)) {
    return(character(0))
  }
  entries <- list.files(patient_dir, full.names = TRUE, recursive = FALSE, all.files = FALSE, no.. = TRUE)
  entries[file.info(entries)$isdir %in% TRUE]
}

find_sv_file <- function(sample_dir, patient_id, sample_id, sv_file_pattern) {
  patt <- fill_placeholders(sv_file_pattern, patient_id, sample_id)
  matched <- find_matching_entries(sample_dir, patt, want_dir = FALSE)
  if (length(matched) == 0) NA_character_ else matched[[1]]
}

find_cn_file <- function(cn_root, patient_id, sample_id, cn_dir_pattern, cn_file_pattern) {
  dir_patt <- fill_placeholders(cn_dir_pattern, patient_id, sample_id)
  matched_dirs <- find_matching_entries(cn_root, dir_patt, want_dir = TRUE)

  if (length(matched_dirs) == 0) {
    return(NA_character_)
  }

  file_patt <- fill_placeholders(cn_file_pattern, patient_id, sample_id)

  for (d in matched_dirs) {
    matched_files <- find_matching_entries(d, file_patt, want_dir = FALSE)
    if (length(matched_files) > 0) {
      return(matched_files[[1]])
    }
  }

  NA_character_
}

build_matrix <- function(patients, sv_root, cn_root, sv_file_pattern, cn_dir_pattern, cn_file_pattern) {
  rows <- vector("list", length = 0)

  for (pid in patients) {
    patient_dir <- file.path(sv_root, pid)
    sample_dirs <- list_sample_dirs(patient_dir)

    if (length(sample_dirs) == 0) {
      next
    }

    for (sdir in sample_dirs) {
      sid <- basename(sdir)
      sv_file <- find_sv_file(sdir, pid, sid, sv_file_pattern)
      if (is.na(sv_file)) {
        next
      }

      cn_file <- find_cn_file(cn_root, pid, sid, cn_dir_pattern, cn_file_pattern)
      rows[[length(rows) + 1]] <- data.frame(
        patientID = pid,
        sampleID = sid,
        SV_file_dir = sv_file,
        CN_file_dir = cn_file,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) {
    data.frame(
      patientID = character(0),
      sampleID = character(0),
      SV_file_dir = character(0),
      CN_file_dir = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    do.call(rbind, rows)
  }
}

main <- function() {
  args <- parse_args()

  if (isTRUE(args$help) || length(args) == 0) {
    print_usage()
    return(invisible(NULL))
  }

  patients_file <- args$patients %||% stop("Missing --patients")
  sv_root <- args$`sv-root` %||% stop("Missing --sv-root")
  cn_root <- args$`cn-root` %||% stop("Missing --cn-root")
  out_file <- args$out %||% stop("Missing --out")
  sv_file_pattern <- args$`sv-file-pattern` %||% "{sampleID}\\.somatic\\.sv\\.bedpe"
  cn_dir_pattern <- args$`cn-dir-pattern` %||% "P_{patientID}_T_{sampleID}_.*"
  cn_file_pattern <- args$`cn-file-pattern` %||% "{sampleID}_copynumber_1_gapfilled_LogR_new\\.txt"

  if (!dir.exists(sv_root)) stop("SV root directory not found: ", sv_root)
  if (!dir.exists(cn_root)) stop("CN root directory not found: ", cn_root)

  patients <- read_patients(patients_file)
  mat <- build_matrix(patients, sv_root, cn_root, sv_file_pattern, cn_dir_pattern, cn_file_pattern)

  out_dir <- dirname(out_file)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  write.table(mat, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)

  message("Wrote matrix: ", out_file)
  message("Rows: ", nrow(mat))
  message("Rows with CN file: ", sum(!is.na(mat$CN_file_dir) & nzchar(mat$CN_file_dir)))
}

main()
