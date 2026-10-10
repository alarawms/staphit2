#!/usr/bin/env python3
"""Split a recombination-masked whole-genome alignment into BEAST2 inputs.

    beast_alignment.py <masked.full.aln> <snps.fasta> <constant_sites.txt> [--drop Reference]

snps.fasta          columns where the samples carry >= 2 different bases (A/C/G/T)
constant_sites.txt  "A C G T" counts of columns where every sample has the same base;
                    BEAST2 needs them (constantSiteWeights) to correct for analysing
                    variable sites only, otherwise rates are overestimated

Columns with N, gaps or masked (recombinant) bases in any sample are left out of both.
Sequences named with --drop (snippy-core's Reference row) are removed first.
"""
import argparse

BASES = 'ACGT'


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.upper())
    return {k: ''.join(v) for k, v in seqs.items()}


def split(seqs):
    names = list(seqs)
    length = len(seqs[names[0]])
    keep, constant = [], dict.fromkeys(BASES, 0)
    for i in range(length):
        column = {seqs[n][i] for n in names}
        if not column <= set(BASES):
            continue
        if len(column) == 1:
            constant[column.pop()] += 1
        else:
            keep.append(i)
    snps = {n: ''.join(seqs[n][i] for i in keep) for n in names}
    return snps, constant


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('alignment')
    ap.add_argument('snps')
    ap.add_argument('constant')
    ap.add_argument('--drop', action='append', default=[])
    a = ap.parse_args(argv)

    seqs = {k: v for k, v in read_fasta(a.alignment).items() if k not in a.drop}
    snps, constant = split(seqs)
    with open(a.snps, 'w') as fh:
        for name, seq in snps.items():
            fh.write(f'>{name}\n{seq}\n')
    with open(a.constant, 'w') as fh:
        fh.write(' '.join(str(constant[b]) for b in BASES) + '\n')
    n_snps = len(next(iter(snps.values()))) if snps else 0
    print(f'{len(snps)} sequences, {n_snps} variable sites, constant A C G T = '
          + ' '.join(str(constant[b]) for b in BASES))


if __name__ == '__main__':
    main()
