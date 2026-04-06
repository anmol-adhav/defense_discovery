# Defense Island Novel Protein Discovery Pipeline

A computational pipeline to identify candidate novel bacterial defense proteins by
mining genomic neighborhoods of known defense systems, filtering prophage contamination,
and prioritizing candidates by multi-system co-occurrence and structural novelty.

## Overview

```
Bacterial genomes (.fna + .gff)
        │
        ▼
PADLOC + DefenseFinder annotation
        │
        ▼
Defense island construction
        │
        ▼
Neighbor protein extraction (±3 genes of known defense proteins)
        │
        ▼
PhiSpy prophage masking
        │
        ▼
MMseqs2 clustering (30% identity, 80% coverage)
        │
        ▼
Multi-system co-occurrence scoring
        │
        ▼
Pfam annotation (hmmscan)
        │
        ▼
Structural prediction (ESMfold) + Foldseek novelty search
```

## Key Results

Two structurally novel, uncharacterized proteins identified from 72 bacterial genomes:

| Candidate | Genomes | System Types | Top System | Pfam | Foldseek |
|-----------|---------|-------------|------------|------|----------|
| `RRWJ01000050.1_27` | 14/72 | 6 | PD-T7-1or5 | none | no hits |
| `RRWJ01000068.1_11` | 19/72 | 4 | AbiE | none | pending |

## Repository Structure

```
defense_discovery/
├── README.md
├── LICENSE
├── .gitignore
├── envs/
│   └── defense_discovery.yml       # conda environment
└── scripts/
    ├── 01_run_pyrodigal.sh         # gene prediction
    ├── 02_run_padloc.sh            # defense system annotation
    ├── 03_run_defensefinder.sh     # defense system annotation
    ├── 04_merge_results.py         # merge PADLOC + DefenseFinder outputs
    ├── 05_build_islands.py         # cluster defense genes into islands
    ├── 06_find_neighbors.py        # extract ±3 gene neighbors
    ├── 07_run_phispy.py            # prophage detection on all genomes
    ├── 08_parse_phispy_logs.py     # extract prophage coordinates
    ├── 09_filter_prophage_neighbors.py
    ├── 10_cluster_neighbors.sh     # MMseqs2 clustering
    ├── 11_summarize_clusters.py    # multi-system scoring
    ├── 12_annotate_candidates.py   # Pfam annotation
    └── 13_extract_top_sequences.py # extract FASTA for structural prediction
```

## Dependencies

- Python ≥ 3.11
- [PADLOC](https://github.com/padlocbio/padloc)
- [DefenseFinder](https://github.com/mdmparis/defense-finder)
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [HMMER 3.4](http://hmmer.org/) + Pfam-A database
- [PhiSpy v5](https://github.com/linsalrob/PhiSpy)
- [Pyrodigal](https://github.com/althonos/pyrodigal)
- Biopython, pandas, numpy, scikit-learn

See `envs/defense_discovery.yml` for the full conda environment.

## Installation

```bash
git clone https://github.com/YOUR_USERNAME/defense_discovery.git
cd defense_discovery
conda env create -f envs/defense_discovery.yml
conda activate defense_discovery

# Download Pfam-A HMM database (~500 MB)
mkdir -p ~/pfam
wget -q https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -P ~/pfam/
gunzip ~/pfam/Pfam-A.hmm.gz
hmmpress ~/pfam/Pfam-A.hmm
```

## Input Data Layout

```
project/
├── accessions.txt          # one NCBI assembly accession per line
├── genomes/                # .fna genome FASTA files
├── gff/                    # .gff gene predictions (from Pyrodigal)
├── proteins/               # .faa protein sequences (from Pyrodigal)
├── gbk/                    # GenBank files (downloaded via ncbi-datasets, for PhiSpy)
├── padloc_out/             # PADLOC output (one subdir per genome)
└── defensefinder_out/      # DefenseFinder output (one subdir per genome)
```

## Running the Pipeline

```bash
# Step 1-3: Gene prediction + defense annotation
bash scripts/01_run_pyrodigal.sh
bash scripts/02_run_padloc.sh
bash scripts/03_run_defensefinder.sh

# Step 4-5: Merge and build islands
python scripts/04_merge_results.py
python scripts/05_build_islands.py

# Step 6: Extract neighbors
python scripts/06_find_neighbors.py

# Step 7-9: Prophage masking
python scripts/07_run_phispy.py
python scripts/08_parse_phispy_logs.py
python scripts/09_filter_prophage_neighbors.py

# Step 10: Cluster neighbors
bash scripts/10_cluster_neighbors.sh

# Step 11-13: Score, annotate, extract
python scripts/11_summarize_clusters.py
python scripts/12_annotate_candidates.py
python scripts/13_extract_top_sequences.py
```

## Output Files

| File | Description |
|------|-------------|
| `all_padloc.csv` | Merged PADLOC hits across all genomes |
| `all_defensefinder.csv` | Merged DefenseFinder hits |
| `defense_islands.csv` | Defense gene clusters per genome |
| `defense_neighbors.csv` | Raw neighbor proteins (±3 genes) |
| `prophage_regions.csv` | PhiSpy prophage coordinates |
| `defense_neighbors_clean.csv` | Neighbors after prophage filtering |
| `candidate_novel_defense_clean.csv` | Ranked candidates by multi-system score |
| `candidate_novel_defense_annotated.csv` | Final candidates with Pfam annotations |
| `top_reps.faa` | Representative sequences for structural prediction |

## Scoring Logic

Candidates are ranked by `n_system_types`: the number of distinct defense system
types (e.g. RM_type_I, CBASS, AbiE) each protein cluster co-occurs with across
all genomes. A high `n_system_types` score with no known Pfam domain is the
primary indicator of a novel defense-associated protein.

Filters applied:
- Present in ≥ 10 genomes (`n_genomes`)
- Co-occurs with ≥ 2 distinct defense system types (`n_system_types`)
- Not located inside a PhiSpy-predicted prophage region
- Not annotated by Pfam (final novelty filter)

## Citation

If you use this pipeline, please cite:
- [PADLOC](https://doi.org/10.1093/nar/gkac400) — Payne et al., 2022
- [DefenseFinder](https://doi.org/10.1038/s41467-022-30269-9) — Tesson et al., 2022
- [MMseqs2](https://doi.org/10.1038/nbt.3988) — Steinegger & Söding, 2017
- [PhiSpy](https://doi.org/10.1093/nar/gks406) — Akhter et al., 2012
- [HMMER](http://hmmer.org) — Eddy, 2011
- [Pyrodigal](https://doi.org/10.21105/joss.04296) — Larralde, 2022

## License

MIT License — see `LICENSE`.
