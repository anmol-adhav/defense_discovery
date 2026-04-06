#!/usr/bin/env bash
# DefenseFinder defense system annotation
set -euo pipefail

PROTEINS_DIR="proteins"
OUT_DIR="defensefinder_out"
THREADS=4

mkdir -p "$OUT_DIR"

for faa in "$PROTEINS_DIR"/*.faa; do
    genome=$(basename "$faa" .faa)
    if [[ -d "$OUT_DIR/$genome" ]]; then
        echo "Skip $genome (exists)"
        continue
    fi
    echo "Running DefenseFinder on $genome..."
    mkdir -p "$OUT_DIR/$genome"
    defense-finder run "$faa" \
        --out-dir "$OUT_DIR/$genome" \
        --workers "$THREADS"
done
echo "DefenseFinder complete."
