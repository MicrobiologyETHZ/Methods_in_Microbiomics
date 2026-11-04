#!/usr/bin/env bash
# Reference solution, Part C: quantifying the catalog in each sample.
# Usage: bash solutions/part_c.sh [threads]
set -euo pipefail
THREADS=${1:-4}
cd "$(dirname "$0")/.."
SAMPLES="metag1 metag2 metag3"

# ---------------------------------------------------------------------------
# Part C: quantification
# ---------------------------------------------------------------------------
mkdir -p alignments
bwa index catalog/catalog.fna 2> /dev/null
for s in $SAMPLES; do
  for r in 1 2; do
    bwa mem -a -t "$THREADS" catalog/catalog.fna data/reads/$s.$r.fq.gz 2> /dev/null \
      | samtools view -F 4 -b - > alignments/$s.r$r.bam
    # read name, flag, gene, CIGAR, NM (edit distance)
    samtools view alignments/$s.r$r.bam \
      | awk -v OFS='\t' '{ nm = "NA"; for (i = 12; i <= NF; i++) if ($i ~ /^NM:i:/) nm = substr($i, 6);
                           print $1, $2, $3, $6, nm }' \
      > alignments/$s.r$r.tsv
  done
done

Rscript solutions/count_genes.R alignments data/reads counts
