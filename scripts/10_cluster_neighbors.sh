#!/usr/bin/env bash
# Cluster neighbor proteins with MMseqs2 at 30% identity / 80% coverage
set -euo pipefail

mmseqs easy-cluster neighbor_proteins_clean.faa neighbor_clusters_clean tmp_clean \
    --min-seq-id 0.30 \
    -c 0.80 \
    --cov-mode 0 \
    --cluster-mode 2 \
    --threads 4 \
    -v 1

echo "Clustering complete."
echo "Clusters: $(cut -f1 neighbor_clusters_clean_cluster.tsv | sort -u | wc -l)"
