#!/usr/bin/env bash
# Gene prediction with Pyrodigal on all genomes
set -euo pipefail

GENOMES_DIR="genomes"
GFF_DIR="gff"
PROTEINS_DIR="proteins"
THREADS=4

mkdir -p "$GFF_DIR" "$PROTEINS_DIR"

for fna in "$GENOMES_DIR"/*.fna; do
    genome=$(basename "$fna" .fna)
    if [[ -f "$PROTEINS_DIR/${genome}.faa" ]]; then
        echo "Skip $genome (exists)"
        continue
    fi
    echo "Processing $genome..."
    pyrodigal -i "$fna" \
              -a "$PROTEINS_DIR/${genome}.faa" \
              -f gff \
              -o "$GFF_DIR/${genome}.gff" \
              -p meta \
              --jobs "$THREADS"
done
echo "Pyrodigal complete."
