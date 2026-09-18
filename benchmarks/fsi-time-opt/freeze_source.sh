#!/bin/bash
# Lightweight source-only archival; never archive HEAD instead of the dirty tree.
set -euo pipefail
variant=${1:?variant name required}
root=${2:?evidence root required}
test ! -e "$root/$variant-source.tar"
mkdir -p "$root"
git ls-files --cached --others --exclude-standard | LC_ALL=C sort -u |
    awk '/^(include\/|solvers\/cpu\/(include\/|src\/|tests\/|Makefile$)|benchmarks\/fsi-time-opt\/|examples\/vascular_flow\/bifurcation_fsi\/)/ && (/\.(hpp|h|cpp|c|sh|py|sbatch|json|mk)$/ || /\/Makefile$/) { if ($0 !~ /benchmarks\/fsi-time-opt\/[0-9]/ && $0 !~ /\/results\//) print }' > "$root/$variant-files.txt"
test -s "$root/$variant-files.txt"
tar -cf "$root/$variant-source.tar" -T "$root/$variant-files.txt"
xargs -d '\n' sha256sum < "$root/$variant-files.txt" > "$root/$variant-source-files.sha256"
sha256sum "$root/$variant-source.tar" > "$root/$variant-source.sha256"
git diff --binary > "$root/$variant-worktree.patch"
git status --short > "$root/$variant-git-status.txt"
git rev-parse HEAD > "$root/$variant-head.txt"
