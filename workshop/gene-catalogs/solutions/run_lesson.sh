#!/usr/bin/env bash
# Reference solution for the gene-catalog lesson: runs Parts A–C in order.
# Usage: bash solutions/run_lesson.sh [threads]
set -euo pipefail
THREADS=${1:-4}
DIR="$(dirname "$0")"
bash "$DIR/part_a.sh" "$THREADS"
bash "$DIR/part_b.sh" "$THREADS"
bash "$DIR/part_c.sh" "$THREADS"
echo "Done. Results in genes/, catalog/, alignments/ and counts/."
