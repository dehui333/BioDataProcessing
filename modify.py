#!/usr/bin/env python

import argparse
from Bio import SeqIO
from functools import partial
from pathlib import Path
import pysam
'''
Various ways to modify the records in a sequence file.
* I think just printing output to stdout will be faster?
* should try to use low level parser from SeqIO such as from Bio.SeqIO.FastaIO import SimpleFastaParser
'''

'''
Rename each sequence in a orderly fashion.
'''    
def ordered_names(input_path, output_path):
    file_type = Path(input_path).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'
    output_name = Path(output_path).stem
    records =  list(SeqIO.parse(input_path, file_type))
    index = 0
    for record in records:
        record.id = output_name + '_' + str(index)
        record.name = ''
        record.description = ''
        index += 1
    SeqIO.write(records, output_path, file_type)

'''
Append a suffix to all sequence names.
'''    
def suffix_names(input_path, output_path, suffix):
    file_type = Path(input_path).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'
    records =  list(SeqIO.parse(input_path, file_type))
    for record in records:
        if record.name == record.id:
            record.name = ''
        if record.description == record.id:
            record.description = ''
        record.id = record.id + suffix
    if file_type in ['fq']:
        file_type = 'fastq'
    SeqIO.write(records, output_path, file_type)

def bam_pair(input_path, output_path):
    samfile = pysam.AlignmentFile(input_path, "r")
    output = pysam.AlignmentFile(output_path, "w", template=samfile)
    for read in samfile.fetch():
        if read.is_paired:
            if read.is_read1:
                read.query_name = read.query_name + '/1'
            if read.is_read2:
                read.query_name = read.query_name + '/2'
        output.write(read)
    output.close()
    samfile.close()

def longest_hp(seq):
    current_char = 'X'
    current_count = 0
    max_count = 0
    for c in seq:
        if c == current_char:
            current_count += 1
        else:
            if current_count > max_count:
                max_count = current_count
            current_count = 1
            current_char = c
    if current_count > max_count:
                max_count = current_count 
    return max_count

def has_long_HP(seq_record, long_len):
    return longest_hp(seq_record) >= long_len

def filter_reads(input_path, output_path, predicate):
    file_type = Path(input_path).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'
    records =  list(SeqIO.parse(input_path, file_type))
    new_records = []
    for record in records:
        if predicate(record):
            new_records.append(record)
    if file_type in ['fq']:
        file_type = 'fastq'
    SeqIO.write(new_records, output_path, file_type)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Some utility functionalities.')
    subparsers = parser.add_subparsers(dest='command')

    order_parser = subparsers.add_parser("order", help='modify sequence names in a fastx file to filename_index.')
    order_parser.add_argument('-i', '--input', type=str, help='path to input fasta/q file.')
    order_parser.add_argument('-o', '--output', type=str, help='path to output fasta/q file.')
    
    suffix_parser = subparsers.add_parser("suffix", help='modify sequence names in a fastx file to filename_suffix.')
    suffix_parser.add_argument('-i', '--input', type=str, help='path to input fasta/q file.')
    suffix_parser.add_argument('-o', '--output', type=str, help='path to output fasta/q file.')
    suffix_parser.add_argument('-s', '--suffix', type=str, help='suffix to add to every name.')

    bam_pair_parser = subparsers.add_parser("bam_pair", help='modify query names in a bam/sam file to indicate pairing.')
    bam_pair_parser.add_argument('-i', '--input', type=str, help='path to input file.')
    bam_pair_parser.add_argument('-o', '--output', type=str, help='path to output file.')

    filter_parser = subparsers.add_parser("filter", help='Filter to retain sequences with HP len >= l.')
    filter_parser.add_argument('-i', '--input', type=str, help='path to input file.')
    filter_parser.add_argument('-o', '--output', type=str, help='path to output file.')
    filter_parser.add_argument('-l', '--len', type=int, help='min HP len to retain')


    args = parser.parse_args()
    if args.command == 'order':
        ordered_names(args.input, args.output)
    if args.command == 'suffix':
        suffix_names(args.input, args.output, args.suffix)   
    if args.command == 'bam_pair':
        bam_pair(args.input, args.output)
    if args.command == 'filter':
        f = partial(has_long_HP, long_len=args.len)
        filter_reads(args.input, args.output, f)