#!/usr/bin/env python3
import argparse
from collections import defaultdict
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--gene-matrix', required=True, help='Gene abundance matrix: GeneID + sample columns')
ap.add_argument('--diamond', required=True, help='DIAMOND NR top-hit table from outfmt 6 qseqid ... staxids sscinames')
ap.add_argument('--output', required=True)
ap.add_argument('--unknown-label', default='Unannotated')
args = ap.parse_args()

gene_mat = pd.read_csv(args.gene_matrix, sep='\t')
if 'GeneID' not in gene_mat.columns:
    raise SystemExit('Gene matrix must contain a GeneID column')

ann = pd.read_csv(
    args.diamond,
    sep='\t',
    header=None,
    names=['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'staxids', 'sscinames']
)
ann = ann.drop_duplicates(subset='qseqid', keep='first')

id2tax = {}
for _, row in ann.iterrows():
    tax = str(row['sscinames']).strip()
    if tax == '' or tax.lower() == 'nan':
        tax = args.unknown_label
    id2tax[str(row['qseqid'])] = tax

sample_cols = [c for c in gene_mat.columns if c != 'GeneID']
taxa_to_values = defaultdict(lambda: [0.0] * len(sample_cols))

for _, row in gene_mat.iterrows():
    gene = str(row['GeneID'])
    tax = id2tax.get(gene, args.unknown_label)
    for i, col in enumerate(sample_cols):
        taxa_to_values[tax][i] += float(row[col])

out = pd.DataFrame.from_dict(taxa_to_values, orient='index', columns=sample_cols)
out.index.name = 'Taxon_Name'
out = out.sort_index()
out = out.loc[out.sum(axis=1) > 0, :]
out.to_csv(args.output, sep='\t')
