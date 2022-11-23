#!/usr/bin/env python

import argparse
from Bio import SeqIO
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(
        description='Get reads of the same names as those in the target set from input set.')
    parser.add_argument('-i',
                        '--input',
                        type=str,
                        help='path of the input FASTA/Q file.',
                        required=True)
    parser.add_argument('-t',
                        '--target',
                        type=str,
                        help='path of the target reads.',
                        required=True)
    parser.add_argument('-o',
                        '--out',
                        type=str,
                        help='path of the output file.',
                        required=True)
    args = parser.parse_args()
    return args
def get_seq_file_type(seq_file_path):
    file_type = Path(seq_file_path).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'
    return file_type

def get_target_reads(input_index, target_index, output_path:str):
    with open(output_path, 'w') as output_handle:
        for target_id in target_index:
            output_handle.write(f'>{target_id}\n')
            output_handle.write(f'{str(input_index[target_id].seq)}\n')
             
    
def main():
    args = parse_args()
    input_index = SeqIO.index(args.input, get_seq_file_type(args.input))
    target_index = SeqIO.index(args.target, get_seq_file_type(args.target))
    get_target_reads(input_index, target_index, args.out)

if __name__ == '__main__':
    main()