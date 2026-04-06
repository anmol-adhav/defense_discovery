"""
Score neighbor protein clusters by multi-system co-occurrence.
Requires: neighbor_clusters_clean_cluster.tsv, defense_neighbors_clean.csv, all_padloc.csv
Outputs:  candidate_novel_defense_clean.csv
"""
import pandas as pd

TOTAL_GENOMES = 72   # update to match your dataset

clusters = pd.read_csv('neighbor_clusters_clean_cluster.tsv', sep='\t',
                       header=None, names=['rep', 'member'])
neighbors = pd.read_csv('defense_neighbors_clean.csv')
padloc = pd.read_csv('all_padloc.csv')

prot_to_system = padloc.set_index('target.name')['system'].to_dict()
neighbors['defense_system'] = neighbors['defense_protein'].map(prot_to_system)
merged = clusters.merge(neighbors, left_on='member', right_on='neighbor_protein')

cluster_sizes = clusters.groupby('rep').size().reset_index(name='size')
genome_counts = merged.groupby('rep')['genome'].nunique().reset_index(name='n_genomes')
system_diversity = merged.groupby('rep')['defense_system'].nunique().reset_index(name='n_system_types')
top_system = merged.groupby('rep')['defense_system'].agg(
    lambda x: x.value_counts().index[0]).reset_index(name='top_system')

summary = cluster_sizes.merge(genome_counts, on='rep') \
                        .merge(system_diversity, on='rep') \
                        .merge(top_system, on='rep')
summary['genome_prevalence'] = summary['n_genomes'] / TOTAL_GENOMES

candidates = summary[
    (summary['n_genomes'] >= 10) &
    (summary['n_system_types'] >= 2)
].sort_values(['n_system_types', 'n_genomes'], ascending=False)

print(f"Candidates: {len(candidates)}")
print(candidates.head(20).to_string())
candidates.to_csv('candidate_novel_defense_clean.csv', index=False)
