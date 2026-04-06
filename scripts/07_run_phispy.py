"""
Run PhiSpy prophage detection on all GenBank files.
Requires: gbk/*.gbk  (download with: datasets download genome accession --inputfile accessions.txt --include gbff)
Outputs:  phispy_out/*/phispy.log
"""
import subprocess
import glob
import os
from concurrent.futures import ProcessPoolExecutor

os.makedirs('phispy_out', exist_ok=True)

def run_phispy(gbk):
    genome = os.path.basename(gbk).replace('.gbk', '')
    out = f'phispy_out/{genome}'
    if os.path.exists(f'{out}/phispy.log'):
        return f"Skip {genome}"
    os.makedirs(out, exist_ok=True)
    r = subprocess.run(
        ['PhiSpy.py', gbk, '-o', out, '--output_choice', '4', '--quiet'],
        capture_output=True, text=True
    )
    return f"{'OK' if r.returncode == 0 else 'FAIL'}: {genome}"

gbks = glob.glob('gbk/*.gbk')
print(f"Running PhiSpy on {len(gbks)} genomes with 4 parallel workers...")
with ProcessPoolExecutor(max_workers=4) as ex:
    for result in ex.map(run_phispy, gbks):
        print(result)
print("PhiSpy complete.")
