#!/usr/bin/env bash
# Reference solution, Part A: gene calling.
# Usage: bash solutions/part_a.sh [threads]
set -euo pipefail
THREADS=${1:-4}
cd "$(dirname "$0")/.."
SAMPLES="metag1 metag2 metag3"

# ---------------------------------------------------------------------------
# Part A: gene calling
# ---------------------------------------------------------------------------
mkdir -p genes
for s in $SAMPLES; do
  prodigal -p meta -c -q \
    -i data/scaffolds/$s.scaffolds.min500.fasta \
    -a genes/$s.faa -d genes/$s.fna -f gff -o genes/$s.gff
done

