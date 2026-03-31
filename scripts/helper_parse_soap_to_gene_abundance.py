#!/usr/bin/env python3
import argparse
import csv
from collections import Counter

ap = argparse.ArgumentParser()
ap.add_argument('--soap-pe', required=True)
ap.add_argument('--soap-se', required=True)
ap.add_argument('--gene-lengths', required=True)
ap.add_argument('--sample', required=True)
ap.add_argument('--output', required=True)
args = ap.parse_args()

gene_lengths = {}
with open(args.gene_lengths) as f:
    next(f)
    for line in f:
        gene, length = line.rstrip('\n').split('\t')
        gene_lengths[gene] = int(length)

def count_refs(path, counter):
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            ref = parts[7]
            counter[ref] += 1

counts = Counter()
count_refs(args.soap_pe, counts)
count_refs(args.soap_se, counts)

rpk = {}
for gene, count in counts.items():
    length = gene_lengths.get(gene)
    if length is None or length <= 0:
        continue
    rpk[gene] = count / (length / 1000.0)

total = sum(rpk.values())
with open(args.output, 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(['GeneID', args.sample])
    for gene in sorted(gene_lengths):
        value = 0.0 if total == 0 else rpk.get(gene, 0.0) / total
        writer.writerow([gene, value])
