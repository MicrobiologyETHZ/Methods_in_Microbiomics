#!/usr/bin/env bash
# Reference solution, Part B: clustering the genes into a catalog.
# Usage: bash solutions/part_b.sh [threads]
set -euo pipefail
THREADS=${1:-4}
cd "$(dirname "$0")/.."
SAMPLES="metag1 metag2 metag3"

# ---------------------------------------------------------------------------
# Part B: clustering (commands as in mmseqs_clustering.rst)
# ---------------------------------------------------------------------------
mkdir -p catalog
cat genes/*.fna | gzip > catalog/genes.fna.gz
cat genes/*.faa | gzip > catalog/genes.faa.gz

cd catalog
mmseqs createdb genes.fna.gz genes.db --dbtype 2 --shuffle 0 > mmseqs.log
mmseqs cluster genes.db genes.cluster tmp \
    --kmer-per-seq-scale 0 \
    --kmer-per-seq 1000 \
    -s 4 \
    --max-seq-len 80000 \
    --remove-tmp-files 0 \
    --cluster-mode 2 \
    --min-seq-id 0.95 \
    --threads "$THREADS" \
    --cov-mode 1 \
    -c 0.9 \
    --spaced-kmer-mode 0 \
    --alignment-mode 3 \
    --cluster-reassign 1 >> mmseqs.log
mmseqs createtsv genes.db genes.db genes.cluster genes.cluster.tsv >> mmseqs.log
mmseqs version > mmseqs.version
cd ..

# Name clusters by descending representative length and write the catalog FASTA files
Rscript solutions/name_clusters.R catalog/genes.cluster.tsv catalog/genes.fna.gz catalog/clusters.tsv
cut -f 2 catalog/clusters.tsv | tail -n +2 > catalog/representatives.txt
seqtk subseq catalog/genes.fna.gz catalog/representatives.txt > catalog/representatives.fna
seqtk subseq catalog/genes.faa.gz catalog/representatives.txt > catalog/representatives.faa
# rename representatives to their cluster names (catalog_000000001, ...)
for ext in fna faa; do
  awk 'NR==FNR { if (FNR > 1) name[$2] = $1; next }
       /^>/ { id = substr($1, 2); print ">" name[id]; next } { print }' \
      catalog/clusters.tsv catalog/representatives.$ext > catalog/catalog.$ext
done

