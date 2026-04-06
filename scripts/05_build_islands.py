"""
Cluster neighboring defense genes into genomic islands.
Requires: all_padloc.csv, gff/*.gff
Outputs:  defense_islands.csv
"""
import pandas as pd
import glob
import os

padloc = pd.read_csv('all_padloc.csv')
known_proteins = set(padloc['target.name'].dropna())

# Load gene positions from GFF
gene_positions = {}
for gff in glob.glob('gff/*.gff'):
    genome = os.path.basename(gff).replace('.gff', '')
    with open(gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            attrs = dict(x.split('=') for x in parts[8].split(';') if '=' in x)
            pid = attrs.get('ID', '')
            gene_positions[pid] = {
                'genome': genome, 'contig': parts[0],
                'start': int(parts[3]), 'end': int(parts[4])
            }

# Build islands: defense genes within 10 kb of each other on same contig
defense_gene_info = []
for pid in known_proteins:
    if pid in gene_positions:
        row = gene_positions[pid].copy()
        row['protein'] = pid
        defense_gene_info.append(row)

df = pd.DataFrame(defense_gene_info).sort_values(['genome', 'contig', 'start'])

islands = []
current_island = []
for _, row in df.iterrows():
    if not current_island:
        current_island = [row]
    else:
        prev = current_island[-1]
        if (row['genome'] == prev['genome'] and
                row['contig'] == prev['contig'] and
                row['start'] - prev['end'] <= 10000):
            current_island.append(row)
        else:
            islands.append(current_island)
            current_island = [row]
if current_island:
    islands.append(current_island)

island_rows = []
for i, island in enumerate(islands):
    for gene in island:
        island_rows.append({
            'island_id': i,
            'genome': gene['genome'],
            'contig': gene['contig'],
            'protein': gene['protein'],
            'start': gene['start'],
            'end': gene['end']
        })

result = pd.DataFrame(island_rows)
result.to_csv('defense_islands.csv', index=False)
print(f"Defense islands: {result['island_id'].nunique()}")
print(f"Defense genes in islands: {len(result)}")
