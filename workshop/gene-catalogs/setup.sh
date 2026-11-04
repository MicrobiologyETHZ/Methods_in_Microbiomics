#!/usr/bin/env bash
# Unpack the workshop reads (run automatically when the Codespace is created).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data/reads
tar --no-same-owner -xzf ../../docs/downloads/metag.tar.gz -C data/reads \
    --strip-components=2 --wildcards 'metag_test/reads/*.fq.gz'
echo "Workshop data ready: $(ls data/reads | wc -l) read files, $(ls data/scaffolds | wc -l) assemblies."
