"""
Extract prophage coordinates from PhiSpy log files.
Requires: phispy_out/*/phispy.log
Outputs:  prophage_regions.csv
"""
import re
import glob
import pandas as pd

records = []
for log in glob.glob('phispy_out/*/phispy.log'):
    genome = log.split('/')[1]
    with open(log) as f:
        for line in f:
            m = re.search(r'INFO\s+(\S+)\s+(\d+)\s+(\d+)\s+\d+\s+Kept', line)
            if m:
                records.append({
                    'genome': genome,
                    'contig': m.group(1),
                    'start': int(m.group(2)),
                    'end': int(m.group(3))
                })

df = pd.DataFrame(records)
print(f"Prophage regions: {len(df)} across {df['genome'].nunique()} genomes")
df.to_csv('prophage_regions.csv', index=False)
