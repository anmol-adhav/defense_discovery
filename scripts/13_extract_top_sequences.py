"""
Extract representative FASTA sequences for top candidate clusters.
Requires: candidate_novel_defense_clean.csv (or _annotated.csv),
          neighbor_clusters_clean_cluster.tsv, proteins/*.faa
Outputs:  top_reps.faa        (one seq per cluster rep — for hmmscan/ESMfold)
          top_candidates.faa  (all members of top clusters)
"""
import pandas as pd
from Bio import SeqIO
import glob

TOP_N = 20   # number of top candidates to extract

try:
    candidates = pd.read_csv('candidate_novel_defense_annotated.csv')
    # Prefer novel candidates first
    candidates = pd.concat([
        candidates[candidates['novel'] == True],
        candidates[candidates['novel'] == False]
    ]).head(TOP_N)
except FileNotFoundError:
    candidates = pd.read_csv('candidate_novel_defense_clean.csv').head(TOP_N)

top_reps = set(candidates['rep'])
clusters = pd.read_csv('neighbor_clusters_clean_cluster.tsv', sep='\t',
                       header=None, names=['rep', 'member'])
top_members = set(clusters[clusters['rep'].isin(top_reps)]['member'])

rep_records, all_records = [], []
for faa in glob.glob('proteins/*.faa'):
    for rec in SeqIO.parse(faa, 'fasta'):
        if rec.id in top_reps:
            rep_records.append(rec)
        if rec.id in top_members:
            all_records.append(rec)

SeqIO.write(rep_records, 'top_reps.faa', 'fasta')
SeqIO.write(all_records, 'top_candidates.faa', 'fasta')
print(f"Rep sequences written:     {len(rep_records)}")
print(f"All member seqs written:   {len(all_records)}")
print()
print("Next steps:")
print("  hmmscan --tblout top_candidates_pfam.txt --cpu 4 --noali ~/pfam/Pfam-A.hmm top_reps.faa")
print("  Submit top_reps.faa to ESMfold: https://esmatlas.com/resources?action=fold")
