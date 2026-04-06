"""
Merge PADLOC and DefenseFinder results across all genomes into single CSVs.
Outputs: all_padloc.csv, all_defensefinder.csv
"""
import pandas as pd
import glob
import os

# Merge PADLOC
padloc_frames = []
for f in glob.glob('padloc_out/*/*_padloc.csv'):
    genome = f.split('/')[1]
    df = pd.read_csv(f)
    df['genome'] = genome
    padloc_frames.append(df)

all_padloc = pd.concat(padloc_frames, ignore_index=True)
all_padloc.to_csv('all_padloc.csv', index=False)
print(f"PADLOC: {len(all_padloc)} hits across {all_padloc['genome'].nunique()} genomes")

# Merge DefenseFinder
df_frames = []
for f in glob.glob('defensefinder_out/*/*.tsv'):
    if 'defense_finder_systems' not in f:
        continue
    genome = f.split('/')[1]
    df = pd.read_csv(f, sep='\t')
    df['genome'] = genome
    df_frames.append(df)

all_df = pd.concat(df_frames, ignore_index=True)
all_df.to_csv('all_defensefinder.csv', index=False)
print(f"DefenseFinder: {len(all_df)} hits across {all_df['genome'].nunique()} genomes")
