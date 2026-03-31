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
ap.add_argument('--fasta', required=True)
ap.add_argument('--output', required=True)
args = ap.parse_args()

with open(args.output, 'w') as out:
    out.write('GeneID\tLength\n')
    for name, seq in read_fasta(args.fasta):
        out.write(f'{name}\t{len(seq)}\n')
