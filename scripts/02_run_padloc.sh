#!/usr/bin/env bash
# PADLOC defense system annotation
set -euo pipefail

PROTEINS_DIR="proteins"
OUT_DIR="padloc_out"
THREADS=4

mkdir -p "$OUT_DIR"

for faa in "$PROTEINS_DIR"/*.faa; do
    genome=$(basename "$faa" .faa)
    if [[ -f "$OUT_DIR/${genome}/${genome}_padloc.csv" ]]; then
        echo "Skip $genome (exists)"
        continue
    fi
    echo "Running PADLOC on $genome..."
    padloc --faa "$faa" --outdir "$OUT_DIR/$genome" --threads "$THREADS"
done
echo "PADLOC complete."
