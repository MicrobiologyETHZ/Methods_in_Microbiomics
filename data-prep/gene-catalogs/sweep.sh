#!/usr/bin/env bash
# Catalog size for a grid of identity and coverage thresholds, for the slider on the clustering page.
# All other flags are as in mmseqs_clustering.rst.
# Usage: bash sweep.sh genes.fna.gz sweep.csv [threads]
set -euo pipefail

GENES=$(realpath "$1")
OUT=$(realpath -m "$2")
THREADS=${3:-4}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

mmseqs createdb "$GENES" genes.db --dbtype 2 --shuffle 0 > sweep.log
echo "min_seq_id,coverage,clusters,genes" > "$OUT"
for id in 0.80 0.85 0.90 0.95 0.97 0.99 1.00; do
  for cov in 0.50 0.70 0.80 0.90 0.95; do
    rm -rf tmp genes.cluster*
    mmseqs cluster genes.db genes.cluster tmp \
        --kmer-per-seq-scale 0 --kmer-per-seq 1000 -s 4 --max-seq-len 80000 \
        --remove-tmp-files 0 --cluster-mode 2 --min-seq-id "$id" --threads "$THREADS" \
        --cov-mode 1 -c "$cov" --spaced-kmer-mode 0 --alignment-mode 3 --cluster-reassign 1 \
        ${MMSEQS_DEV_ARGS:-} >> sweep.log
    mmseqs createtsv genes.db genes.db genes.cluster genes.cluster.tsv >> sweep.log
    echo "$id,$cov,$(cut -f 1 genes.cluster.tsv | sort -u | wc -l),$(wc -l < genes.cluster.tsv)" >> "$OUT"
  done
done
