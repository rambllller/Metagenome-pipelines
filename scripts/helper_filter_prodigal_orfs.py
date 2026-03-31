#!/usr/bin/env python3
import argparse
from pathlib import Path

def read_fasta(path):
    name = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name is not None:
            yield name, ''.join(seq)

def write_fasta(records, path):
    with open(path, 'w') as out:
        for name, seq in records:
            out.write(f'>{name}\n')
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + '\n')

ap = argparse.ArgumentParser()
ap.add_argument('--nucleotide', required=True)
ap.add_argument('--protein', required=True)
ap.add_argument('--min-nt-len', type=int, default=100)
ap.add_argument('--out-nucleotide', required=True)
ap.add_argument('--out-protein', required=True)
args = ap.parse_args()

keep = set()
nt_records = []
for name, seq in read_fasta(args.nucleotide):
    if len(seq) >= args.min_nt_len:
        keep.add(name)
        nt_records.append((name, seq))
write_fasta(nt_records, args.out_nucleotide)

aa_records = [(name, seq) for name, seq in read_fasta(args.protein) if name in keep]
write_fasta(aa_records, args.out_protein)
