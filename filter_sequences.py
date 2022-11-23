#!/usr/bin/env python

import argparse
from Bio import SeqIO
from pathlib import Path

def get_seq_file_type(seq_file_path):
    file_type = Path(seq_file_path).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'
    return file_type
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter away sequences shorter than the required length.')
    parser.add_argument('-i', type=str, help='path of the fastx file.')
    parser.add_argument('-l', '--len', type=int, help='Required minimum length.')
    args = parser.parse_args()
    file_type = get_seq_file_type(args.i)
    for record in SeqIO.parse(args.i, file_type):
        if len(record) >= args.len:
            print('>' + record.id)
            print(str(record.seq))