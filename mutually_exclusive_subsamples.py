#!/usr/bin/env python

from Bio import SeqIO
import argparse
from pathlib import Path
from random import shuffle
def parse_args():
    parser = argparse.ArgumentParser(description='Take two mutually exclusive random subsamples from a set of sequences.')
    parser.add_argument('-i', type=str, required=True, help='path of the input fasta/q file.')
    parser.add_argument('-n1', type=int, required=True, help='number to sample for the first subsample.')
    parser.add_argument('-n2', type=int, required=True, help='number to sample for the second subsample.')
    parser.add_argument('-o1', type=str, required=True, help='path of the first output fasta/q file.')
    parser.add_argument('-o2', type=str, required=True, help='path of the second output fasta/q file.')
   
    args = parser.parse_args()
    return args

def get_seq_file_type(seq_file_path):
    file_type = Path(seq_file_path).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'
    return file_type


def main():
    args = parse_args()
    file_type = get_seq_file_type(args.i)
    records = list(SeqIO.parse(args.i, file_type))
    shuffle(records)
    with open(args.o1, 'w') as out1:
        SeqIO.write(records[:args.n1], out1, file_type)
    with open(args.o2, 'w') as out2:
        SeqIO.write(records[args.n1:args.n1+args.n2], out2, file_type)
    
if __name__ == '__main__':
    main()
    