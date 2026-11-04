#!/usr/bin/env bash
# Build all data for the gene-catalog lesson from the raw downloads in docs/downloads/.
#
#   workshop/gene-catalogs/data/   inputs learners use in the Codespace (Parts A–C)
#   site/lessons/gene-catalogs/data/  small precomputed tables for the web pages (Parts A–D)
#
# Requires: metaspades.py, prodigal, mmseqs (Release 18-8cc5c), bwa, samtools, seqtk, R.
# Usage: bash data-prep/gene-catalogs/prepare.sh [threads]
set -euo pipefail

THREADS=${1:-4}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
WORKSHOP=$REPO/workshop/gene-catalogs
SITE_DATA=$REPO/site/lessons/gene-catalogs/data
BUILD=$(mktemp -d)
trap 'rm -rf "$BUILD"' EXIT
SAMPLES="metag1 metag2 metag3"

mkdir -p "$WORKSHOP/data/reads" "$WORKSHOP/data/scaffolds" "$SITE_DATA"

# 1. Toy metagenomes: reads (already quality controlled) and assemblies
tar --no-same-owner -xzf "$REPO/docs/downloads/metag.tar.gz" -C "$BUILD"
cp "$BUILD"/metag_test/reads/*.fq.gz "$WORKSHOP/data/reads/"
for s in $SAMPLES; do
  metaspades.py -t "$THREADS" -m 12 --only-assembler \
    --pe1-1 "$BUILD/metag_test/reads/$s.1.fq.gz" --pe1-2 "$BUILD/metag_test/reads/$s.2.fq.gz" \
    -o "$BUILD/asm/$s" > "$BUILD/asm_$s.log" 2>&1
  python "$REPO/docs/scripts/scaffold_filter.py" "$s" scaffolds \
    "$BUILD/asm/$s/scaffolds.fasta" "$BUILD/asm/$s" META
  cp "$BUILD/asm/$s/$s.scaffolds.min500.fasta" "$WORKSHOP/data/scaffolds/"
done

# 2. Run the full lesson once to produce reference outputs for the check script
REF=$BUILD/reference
mkdir -p "$REF"
cp -r "$WORKSHOP/data" "$WORKSHOP/solutions" "$REF/"
# MMSEQS_DEV_ARGS is for development on machines with < 16 GB RAM only (e.g. "-k 13").
# The official data must be built without it; provenance.json records which was used.
if [ -n "${MMSEQS_DEV_ARGS:-}" ]; then
  echo "WARNING: building PROVISIONAL data with extra mmseqs args: $MMSEQS_DEV_ARGS" >&2
  sed -i "s|--cluster-reassign 1|--cluster-reassign 1 $MMSEQS_DEV_ARGS|" "$REF/solutions/part_b.sh"
fi
bash "$REF/solutions/run_lesson.sh" "$THREADS"
Rscript "$REPO/data-prep/gene-catalogs/reference_summary.R" "$REF" "$WORKSHOP/reference.tsv"

# 3. Tables for the web pages
Rscript "$REPO/data-prep/gene-catalogs/site_tables.R" "$REF" "$SITE_DATA"
bash "$REPO/data-prep/gene-catalogs/sweep.sh" "$REF/catalog/genes.fna.gz" "$SITE_DATA/sweep.csv" "$THREADS"

# 4. Tara Oceans subset for Part D
tar --no-same-owner -xzf "$REPO/docs/downloads/metat_tutorial.tar.gz" -C "$BUILD"
Rscript "$REPO/data-prep/gene-catalogs/tara_tables.R" "$BUILD/metat_tutorial" "$SITE_DATA"

cat > "$SITE_DATA/provenance.json" <<EOF
{
  "built": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "mmseqs_version": "$(cat "$REF/catalog/mmseqs.version")",
  "mmseqs_dev_args": "${MMSEQS_DEV_ARGS:-}",
  "provisional": $([ -n "${MMSEQS_DEV_ARGS:-}" ] && echo true || echo false)
}
EOF
cp "$SITE_DATA/provenance.json" "$WORKSHOP/provenance.json"

echo "Data written to $WORKSHOP/data, $WORKSHOP/reference.tsv and $SITE_DATA"
