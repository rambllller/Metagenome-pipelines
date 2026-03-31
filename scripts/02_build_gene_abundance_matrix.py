#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--input-dir', required=True)
ap.add_argument('--output', required=True)
args = ap.parse_args()

input_dir = Path(args.input_dir)
files = sorted(input_dir.glob('*.gene_abundance.tsv'))
if not files:
    raise SystemExit(f'No *.gene_abundance.tsv files found in {input_dir}')

merged = None
for fp in files:
    df = pd.read_csv(fp, sep='\t')
    if merged is None:
        merged = df
    else:
        merged = merged.merge(df, on='GeneID', how='outer')

merged = merged.fillna(0)
merged.to_csv(args.output, sep='\t', index=False)
