#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# install_shatterseek_RE.R - install ShatterSeek inside a restricted environment
# (e.g. the GEL research environment) where CRAN/GitHub are not reachable and
# the dependencies live in shared library trees.
#
# 1. fetch the source with the raw-URL downloader (no git / API needed):
#      bash gh_raw_fetch.sh parklab/ShatterSeek -o ShatterSeek
# 2. install:
#      Rscript install_shatterseek_RE.R [--src ShatterSeek] [--lib <personal_lib>] [--extra-libs <dir1,dir2>]
#
# Defaults:
#   --src        ./ShatterSeek   (directory produced by gh_raw_fetch.sh, or a .tar.gz)
#   --lib        first writable entry of .libPaths(), else R_LIBS_USER (created if missing)
#   --extra-libs the shared GEL trees, appended to .libPaths() so dependencies are found:
#                /tools/aws-workspace-ubuntu-apps/ce/R/4.5.3 , /tools/aws-workspace-apps/ce/R/4.2.1/
#                (non-existent directories are skipped silently)
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
opt <- list(src = "ShatterSeek", lib = NULL,
            `extra-libs` = "/tools/aws-workspace-ubuntu-apps/ce/R/4.5.3,/tools/aws-workspace-apps/ce/R/4.2.1/")
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[i])
  if (key == "help") { cat(readLines(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]), n = 20), sep = "\n"); quit(status = 0) }
  if (i == length(args)) stop("missing value for --", key)
  opt[[key]] <- args[i + 1]; i <- i + 2
}

# --- library paths ---------------------------------------------------------
extra <- trimws(strsplit(opt$`extra-libs`, ",")[[1]])
extra <- extra[nzchar(extra) & dir.exists(extra)]
.libPaths(c(.libPaths(), extra))          # same effect as the usual .libPaths(c(.libPaths(), "/tools/...")) lines

lib <- opt$lib
if (is.null(lib)) {
  writable <- .libPaths()[file.access(.libPaths(), 2) == 0]
  lib <- if (length(writable)) writable[1] else Sys.getenv("R_LIBS_USER")
}
if (!nzchar(lib)) lib <- file.path("~", "R", paste0(R.version$platform, "-library"), paste(R.version$major, sub("\\..*", "", R.version$minor), sep = "."))
lib <- path.expand(lib)
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

cat("R", R.version.string, "\n")
cat("install into:", lib, "\n")
cat("library search path:\n"); cat(paste0("  ", .libPaths()), sep = "\n")

# --- dependencies ----------------------------------------------------------
deps <- c("methods", "BiocGenerics", "foreach", "graph", "S4Vectors", "GenomicRanges", "IRanges", "MASS", "ggplot2", "grid", "gridExtra")
have <- vapply(deps, requireNamespace, logical(1), quietly = TRUE)
cat("\ndependencies:\n"); cat(sprintf("  %-14s %s", deps, ifelse(have, "ok", "MISSING")), sep = "\n")
if (!all(have)) stop("missing dependencies: ", paste(deps[!have], collapse = ", "),
                     "\nAdd the library tree that contains them with --extra-libs, or install them first.")

# --- install ---------------------------------------------------------------
src <- normalizePath(opt$src, mustWork = FALSE)
if (!file.exists(src)) stop("source not found: ", src, "\nFetch it with:  bash gh_raw_fetch.sh parklab/ShatterSeek -o ShatterSeek")
if (dir.exists(src) && !file.exists(file.path(src, "DESCRIPTION"))) {
  inner <- list.files(src, pattern = "^DESCRIPTION$", recursive = TRUE, full.names = TRUE)   # e.g. ShatterSeek-master/
  if (length(inner)) src <- dirname(inner[1]) else stop("no DESCRIPTION under ", src)
}
cat("\ninstalling from:", src, "\n")
install.packages(src, repos = NULL, type = "source", lib = lib, dependencies = FALSE)

# --- verify ----------------------------------------------------------------
ok <- suppressWarnings(requireNamespace("ShatterSeek", lib.loc = .libPaths(), quietly = TRUE))
if (!ok) stop("ShatterSeek did not install correctly")
cat("\nShatterSeek", as.character(packageVersion("ShatterSeek")), "installed in", find.package("ShatterSeek"), "\n")
cat("Before running the scripts in a new session use the same .libPaths():\n")
cat(sprintf('  .libPaths(c("%s", .libPaths()%s))\n', lib, if (length(extra)) paste0(', "', paste(extra, collapse = '", "'), '"') else ""))
