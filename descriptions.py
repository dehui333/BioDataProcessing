#!/usr/bin/env python

import argparse
from Bio import SeqIO
from pathlib import Path

'''
Parse a description string in key=value key2=value ... format to dictionary
'''
def parse_description(description_string):
    pairs = description_string.split()
    eq_indices = [x.find('=') for x in pairs]
    return {pairs[i][:eq_indices[i]] : pairs[i][eq_indices[i] + 1:] for i in range(len(pairs))}
'''
Convert a description dictionary to a string
'''
def convert_to_string(description_dict):
    output = ''
    for key, value in description_dict.items():
        output += key + '=' + value + ' '
    return output[:-1]

def append_KV_pair_to_record(record, key, value):
    record.description =  key + '=' + value  + record.description[len(record.id):]
    return record

'''
Add info to description of fasta/q sequences in the format of key=value.
'''

def append_KV_pair(seq_iter, key, value):
    for record in seq_iter:
        record.description =  key + '=' + value  + record.description[len(record.id):]
        yield record

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Append a key value pair to the description of fasta/q sequences.')
    parser.add_argument('-i', '--input', type=str, help='Path to input reads.')
    parser.add_argument('-o', '--output', type=str, help='Path of the output file.')
    parser.add_argument('-k', '--key', type=str, help='The key.')
    parser.add_argument('-v', '--value', type=str, help='The value.')
    
    args = parser.parse_args()

    file_type = Path(args.input).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'

    SeqIO.write(append_KV_pair(SeqIO.parse(args.input, file_type), args.key, args.value), args.output, file_type)
    