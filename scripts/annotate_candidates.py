"""
Add Pfam best-hit annotation to candidate clusters.
Requires: candidate_novel_defense_clean.csv, top20_pfam.txt (hmmscan tblout)
Outputs:  candidate_novel_defense_annotated.csv
"""
import pandas as pd

candidates = pd.read_csv('candidate_novel_defense_clean.csv')

pfam_hits = {}
with open('top20_pfam.txt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        query, domain, evalue = parts[2], parts[0], float(parts[4])
        if query not in pfam_hits or evalue < pfam_hits[query][1]:
            pfam_hits[query] = (domain, evalue)

candidates['pfam_hit'] = candidates['rep'].map(
    lambda r: pfam_hits.get(r, ('none', 999))[0])
candidates['pfam_evalue'] = candidates['rep'].map(
    lambda r: pfam_hits.get(r, ('', 999))[1])
candidates['novel'] = candidates['pfam_hit'] == 'none'

print(candidates[['rep', 'size', 'n_genomes', 'n_system_types',
                   'top_system', 'pfam_hit', 'pfam_evalue', 'novel']].to_string())
candidates.to_csv('candidate_novel_defense_annotated.csv', index=False)
