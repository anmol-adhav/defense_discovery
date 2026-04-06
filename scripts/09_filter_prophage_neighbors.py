"""
Remove neighbor proteins located inside PhiSpy-predicted prophage regions.
Requires: prophage_regions.csv, defense_neighbors.csv, gff/*.gff
Outputs:  defense_neighbors_clean.csv, neighbor_proteins_clean.faa
"""
import pandas as pd
import glob
from Bio import SeqIO

prophages = pd.read_csv('prophage_regions.csv')
print(f"Prophage regions: {len(prophages)} across {prophages['genome'].nunique()} genomes")

gene_coords = []
for gff in glob.glob('gff/*.gff'):
    genome = gff.replace('gff/', '').replace('.gff', '')
    with open(gff) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            attrs = dict(x.split('=') for x in parts[8].split(';') if '=' in x)
            gene_coords.append({
                'genome': genome, 'protein': attrs.get('ID', ''),
                'contig': parts[0], 'start': int(parts[3]), 'end': int(parts[4])
            })

genes_df = pd.DataFrame(gene_coords)
merged = genes_df.merge(prophages, on=['genome', 'contig'], how='left', suffixes=('', '_p'))
prophage_proteins = set(merged[
    (merged['start'] >= merged['start_p']) &
    (merged['end'] <= merged['end_p'])
]['protein'])

print(f"Prophage proteins flagged: {len(prophage_proteins)}")

neighbors = pd.read_csv('defense_neighbors.csv')
print(f"Neighbors before filter: {len(neighbors)}")
clean = neighbors[~neighbors['neighbor_protein'].isin(prophage_proteins)]
print(f"Neighbors after filter:  {len(clean)}")
clean.to_csv('defense_neighbors_clean.csv', index=False)

# Extract FASTA sequences for clean neighbors
wanted = set(clean['neighbor_protein'])
records = []
for faa in glob.glob('proteins/*.faa'):
    for rec in SeqIO.parse(faa, 'fasta'):
        if rec.id in wanted:
            records.append(rec)
SeqIO.write(records, 'neighbor_proteins_clean.faa', 'fasta')
print(f"Sequences written: {len(records)}")
