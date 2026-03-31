#!/usr/bin/env python3
import argparse

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

ap = argparse.ArgumentParser()
ap.add_argument('--id-fasta', required=True)
ap.add_argument('--source-fasta', required=True)
ap.add_argument('--output-fasta', required=True)
args = ap.parse_args()

ids = {name for name, _ in read_fasta(args.id_fasta)}
with open(args.output_fasta, 'w') as out:
    for name, seq in read_fasta(args.source_fasta):
        if name in ids:
            out.write(f'>{name}\n')
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + '\n')
