#!/usr/bin/env python

import argparse
from Bio import SeqIO


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter away contigs shorter than the required length.')
    parser.add_argument('-i', '--assm', type=str, help='path of the assembly.')
    parser.add_argument('-l', '--len', type=int, help='Required minimum length.')
    args = parser.parse_args()
    for record in SeqIO.parse(args.assm, 'fasta'):
        if len(record) >= args.len:
            print('>' + record.id)
            print(record.seq)