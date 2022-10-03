#!/usr/bin/env python

import argparse
from Bio import SeqIO
import os
from pathlib import Path

'''
Append a suffix to all sequence ids.
'''    
def add_suffix_to_ids(input_path, output_path, suffix):
    file_type = Path(input_path).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'
    records =  SeqIO.parse(input_path, file_type)
    with open(output_path, 'a') as output_handle: 
        for record in records:
            record.description = record.description[len(record.id)+1:]
            record.name = ''
            record.id = record.id + '_' + suffix
            SeqIO.write(record, output_handle, file_type)
'''
Separate reads into different files based on haplotype, return the file names.
'''
def separate_reads_to_files(reads_path, name1, name2, output_dir):
    file_type = Path(reads_path).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'

    out1 = name1 + '.' + file_type
    out2 = name2 + '.' + file_type
    out1 = output_dir + '/' + Path(reads_path).stem + '-' + out1
    out2 = output_dir + '/' + Path(reads_path).stem + '-' + out2
    if os.path.isfile(out1) and os.path.isfile(out2):
        print(f'{out1} and {out2} already exists. Using existing ones.', file=sys.stderr)
        return out1, out2
    elif os.path.isfile(out1) or os.path.isfile(out2):
        print(f'One of {out1} and {out2} already exists!', file=sys.stderr)
        exit(1)
    with open(out1, 'a') as h1_out, open(out2, 'a') as h2_out:
        for record in SeqIO.parse(reads_path, file_type):
            hap_number = record.id[-1]
            if hap_number == '1':
                SeqIO.write(record, h1_out, file_type)
            elif hap_number == '2':
                SeqIO.write(record, h2_out, file_type)
            else:
                raise ValueError(f'{hap_number} is not an expected haplotype tag!')
    return out1, out2

def main():
    parser = argparse.ArgumentParser(description='Append a suffix to the ids of sequences.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to input reads.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path of the output file.')
    parser.add_argument('-s', '--suffix', type=str, choices = ['1', '2'], required=True, help='The suffix to append.')
    args = parser.parse_args()
    add_suffix_to_ids(args.input, args.output, args.suffix)
if __name__ == '__main__':
    main()