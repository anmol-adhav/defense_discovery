"""
Add Pfam best-hit annotation to candidate clusters and flag unannotated ones.
Requires: candidate_novel_defense_clean.csv, top_candidates_pfam.txt (hmmscan --tblout)
Outputs:  candidate_novel_defense_annotated.csv

Run hmmscan first:
    python 13_extract_top_sequences.py
    hmmscan --tblout top_candidates_pfam.txt --cpu 4 --noali ~/pfam/Pfam-A.hmm top_reps.faa
"""
import pandas as pd
import sys

PFAM_TBLOUT = 'top_candidates_pfam.txt'

candidates = pd.read_csv('candidate_novel_defense_clean.csv')

pfam_hits = {}
try:
    with open(PFAM_TBLOUT) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            query, domain, evalue = parts[2], parts[0], float(parts[4])
            if query not in pfam_hits or evalue < pfam_hits[query][1]:
                pfam_hits[query] = (domain, evalue)
except FileNotFoundError:
    print(f"WARNING: {PFAM_TBLOUT} not found. Run hmmscan first.")
    sys.exit(1)

candidates['pfam_hit']    = candidates['rep'].map(lambda r: pfam_hits.get(r, ('none', 999))[0])
candidates['pfam_evalue'] = candidates['rep'].map(lambda r: pfam_hits.get(r, ('', 999))[1])
candidates['novel']       = candidates['pfam_hit'] == 'none'

novel = candidates[candidates['novel']]
print(f"Total candidates: {len(candidates)}")
print(f"Novel (no Pfam hit): {len(novel)}")
print()
print(candidates[['rep', 'size', 'n_genomes', 'n_system_types',
                   'top_system', 'pfam_hit', 'pfam_evalue', 'novel']].to_string())

candidates.to_csv('candidate_novel_defense_annotated.csv', index=False)
