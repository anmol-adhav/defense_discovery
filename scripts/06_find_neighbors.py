"""
Extract proteins in the ±3 gene neighborhood of known defense proteins.
Requires: all_padloc.csv, gff/*.gff
Outputs:  defense_neighbors.csv
"""
import pandas as pd
import glob
import os

padloc = pd.read_csv('all_padloc.csv')
known_proteins = set(padloc['target.name'].dropna())

neighbors = []
for gff in glob.glob('gff/*.gff'):
    genome = os.path.basename(gff).replace('.gff', '')
    proteins = []
    with open(gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            attrs = dict(x.split('=') for x in parts[8].split(';') if '=' in x)
            pid = attrs.get('ID', '')
            proteins.append((parts[0], int(parts[3]), int(parts[4]), pid))
    proteins.sort(key=lambda x: (x[0], x[1]))

    for i, (seqid, start, end, pid) in enumerate(proteins):
        if pid in known_proteins:
            for j in range(max(0, i - 3), min(len(proteins), i + 4)):
                npid = proteins[j][3]
                if npid not in known_proteins:
                    neighbors.append({
                        'genome': genome,
                        'neighbor_protein': npid,
                        'defense_protein': pid,
                        'distance': j - i
                    })

df = pd.DataFrame(neighbors).drop_duplicates('neighbor_protein')
df.to_csv('defense_neighbors.csv', index=False)
print(f"Unique neighbor proteins: {len(df)}")
