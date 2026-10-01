#!/usr/bin/env bash
# gh_raw_fetch.sh - download files from any PUBLIC GitHub repo using only the
# /raw/ address (https://github.com/<owner>/<repo>/raw/<branch>/<path>), which is
# the one route that works from restricted environments. No git, no API.
#
# Usage:
#   gh_raw_fetch.sh <owner>/<repo>[@<branch>] [-o <outdir>] <path> [<path> ...]
#   gh_raw_fetch.sh <owner>/<repo>[@<branch>] [-o <outdir>] -l <listfile>   # one repo path per line
#   gh_raw_fetch.sh <owner>/<repo>[@<branch>] [-o <outdir>]                 # no paths: fetch FILES.txt from the
#                                                                           # repo root and download everything in it
# Examples:
#   gh_raw_fetch.sh Tim-Yu/chromothripsis-shatterseek scripts/shatterseek_lib.R README.md
#   gh_raw_fetch.sh Tim-Yu/chromothripsis-shatterseek@main -o chromothripsis
#   gh_raw_fetch.sh parklab/ShatterSeek -l shatterseek_files.txt
# Branch defaults to main (falls back to master). Uses curl if present, else wget.
set -euo pipefail

usage() { sed -n '2,16p' "$0"; exit 1; }
[ $# -ge 1 ] || usage

repo_spec=$1; shift
repo=${repo_spec%@*}
branch=main; [[ "$repo_spec" == *@* ]] && branch=${repo_spec#*@}
[[ "$repo" == */* ]] || { echo "repo must be <owner>/<repo>" >&2; usage; }

outdir=.
listfile=""
paths=()
while [ $# -gt 0 ]; do
  case "$1" in
    -o) outdir=$2; shift 2 ;;
    -l) listfile=$2; shift 2 ;;
    -h|--help) usage ;;
    *) paths+=("$1"); shift ;;
  esac
done

if command -v curl >/dev/null 2>&1; then
  fetch() { curl -fsL --retry 3 -o "$2" "$1" 2>/dev/null; }
elif command -v wget >/dev/null 2>&1; then
  fetch() { wget -q -O "$2" "$1" 2>/dev/null; }
else
  echo "need curl or wget" >&2; exit 1
fi

raw_url() { echo "https://github.com/${repo}/raw/${1}/${2}"; }

# branch: try the requested one on a probe file, fall back to master when 'main' was only the default
pick_branch() {
  fetch "$(raw_url "$branch" "$1")" /dev/null && return 0
  if [ "$branch" = main ] && fetch "$(raw_url master "$1")" /dev/null; then branch=master; return 0; fi
  return 1
}

if [ -n "$listfile" ]; then
  mapfile -t paths < <(grep -v '^\s*#' "$listfile" | sed '/^\s*$/d')
  pick_branch "${paths[0]}" || true
elif [ ${#paths[@]} -eq 0 ]; then
  pick_branch FILES.txt || { echo "no paths given and no FILES.txt in ${repo} (main/master)" >&2; exit 1; }
  tmp=$(mktemp); fetch "$(raw_url "$branch" FILES.txt)" "$tmp"
  mapfile -t paths < <(grep -v '^\s*#' "$tmp" | sed '/^\s*$/d'); rm -f "$tmp"
else
  pick_branch "${paths[0]}" || true
fi

ok=0; fail=0
for p in "${paths[@]}"; do
  dest="${outdir}/${p}"
  mkdir -p "$(dirname "$dest")"
  if fetch "$(raw_url "$branch" "$p")" "$dest"; then
    echo "ok    ${repo}@${branch}:${p} -> ${dest}"; ok=$((ok+1))
  else
    echo "FAIL  ${repo}@${branch}:${p}" >&2; rm -f "$dest"; fail=$((fail+1))
  fi
done
echo "downloaded ${ok}, failed ${fail}"
[ "$fail" -eq 0 ]
