#!/usr/bin/env bash
# gh_raw_fetch.sh - download files (or a whole repository) from any PUBLIC GitHub
# repo using only github.com URLs: files come from the /raw/ address
# (https://github.com/<owner>/<repo>/raw/<branch>/<path>) and the file list is
# read from the repository's tree pages (https://github.com/<owner>/<repo>/tree/<branch>/<dir>
# requested as JSON). No git, no api.github.com, no index file needed.
#
# Usage:
#   gh_raw_fetch.sh <owner>/<repo>[@<branch>] [-o <outdir>]                  # whole repository
#   gh_raw_fetch.sh <owner>/<repo>[@<branch>] [-o <outdir>] <path> [<path>...] # selected files
#   gh_raw_fetch.sh <owner>/<repo>[@<branch>] [-o <outdir>] -l <listfile>    # paths listed one per line
#   gh_raw_fetch.sh <owner>/<repo>[@<branch>] [-o <outdir>] -d <dir>         # one sub-directory, recursively
# Examples:
#   gh_raw_fetch.sh Tim-Yu/chromothripsis-shatterseek -o chromothripsis-shatterseek
#   gh_raw_fetch.sh parklab/ShatterSeek@master -d R -o ShatterSeek
#   gh_raw_fetch.sh Tim-Yu/chromothripsis-shatterseek scripts/shatterseek_lib.R README.md
# <branch> may also be a tag or commit SHA; default main, falling back to master.
# Uses curl if present, otherwise wget. Requires only bash, grep, sed.
set -euo pipefail

usage() { sed -n '2,17p' "$0"; exit 1; }
[ $# -ge 1 ] || usage

repo_spec=$1; shift
repo=${repo_spec%@*}
branch=main; [[ "$repo_spec" == *@* ]] && branch=${repo_spec#*@}
[[ "$repo" == */* ]] || { echo "repo must be <owner>/<repo>" >&2; usage; }

outdir=. ; listfile="" ; subdir="" ; paths=()
while [ $# -gt 0 ]; do
  case "$1" in
    -o) outdir=$2; shift 2 ;;
    -l) listfile=$2; shift 2 ;;
    -d) subdir=${2%/}; shift 2 ;;
    -h|--help) usage ;;
    -*) echo "unknown option $1" >&2; usage ;;
    *) paths+=("$1"); shift ;;
  esac
done

if command -v curl >/dev/null 2>&1; then
  fetch()      { curl -fsL --retry 3 -o "$2" "$1" 2>/dev/null; }
  fetch_json() { curl -fsL --retry 3 -H "Accept: application/json" "$1" 2>/dev/null; }
elif command -v wget >/dev/null 2>&1; then
  fetch()      { wget -q -O "$2" "$1" 2>/dev/null; }
  fetch_json() { wget -q -O - --header="Accept: application/json" "$1" 2>/dev/null; }
else
  echo "need curl or wget" >&2; exit 1
fi

raw_url()  { echo "https://github.com/${repo}/raw/${branch}/${1}"; }
tree_url() { echo "https://github.com/${repo}/tree/${branch}/${1}"; }

# List one directory from the tree page JSON: prints "<type>\t<path>" per entry.
list_dir() {
  fetch_json "$(tree_url "$1")" \
    | grep -o '"path":"[^"]*","contentType":"[a-z_]*"' \
    | sed 's/^"path":"\(.*\)","contentType":"\(.*\)"$/\2\t\1/'
}

# Recursive walk; prints file paths.
walk_tree() {
  local dir=$1 type path
  while IFS=$'\t' read -r type path; do
    [ -n "$path" ] || continue
    case "$type" in
      directory) walk_tree "$path" ;;
      file|symlink_file) echo "$path" ;;
      *) echo "skip  ${type}: ${path}" >&2 ;;        # submodules etc.
    esac
  done < <(list_dir "$dir")
}

branch_exists() { fetch_json "$(tree_url "")" | grep -q '"contentType"'; }

# resolve the branch when 'main' was only the default
if ! branch_exists; then
  if [ "$branch" = main ] && { branch=master; branch_exists; }; then :; else
    echo "cannot read tree of ${repo}@${branch} (private repo, wrong branch, or github.com not reachable)" >&2
    [ ${#paths[@]} -gt 0 ] || [ -n "$listfile" ] || exit 1   # explicit paths can still be tried on the given branch
    branch=${repo_spec#*@}; [[ "$repo_spec" == *@* ]] || branch=main
  fi
fi

if [ -n "$listfile" ]; then
  mapfile -t paths < <(grep -v '^\s*#' "$listfile" | sed '/^\s*$/d')
elif [ -n "$subdir" ]; then
  mapfile -t paths < <(walk_tree "$subdir")
elif [ ${#paths[@]} -eq 0 ]; then
  mapfile -t paths < <(walk_tree "")
  if [ ${#paths[@]} -eq 0 ]; then                    # fallback: an index file in the repo root
    tmp=$(mktemp)
    if fetch "$(raw_url FILES.txt)" "$tmp"; then mapfile -t paths < <(grep -v '^\s*#' "$tmp" | sed '/^\s*$/d'); fi
    rm -f "$tmp"
  fi
  [ ${#paths[@]} -gt 0 ] || { echo "could not list files of ${repo}@${branch}" >&2; exit 1; }
fi

echo "${repo}@${branch}: ${#paths[@]} file(s) -> ${outdir}"
ok=0; fail=0
for p in "${paths[@]}"; do
  dest="${outdir}/${p}"
  mkdir -p "$(dirname "$dest")"
  if fetch "$(raw_url "$p")" "$dest"; then
    echo "ok    ${p}"; ok=$((ok+1))
  else
    echo "FAIL  ${p}" >&2; rm -f "$dest"; fail=$((fail+1))
  fi
done
echo "downloaded ${ok}, failed ${fail}"
[ "$fail" -eq 0 ]
