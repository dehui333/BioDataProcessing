#!/usr/bin/env python

import argparse
from Bio import SeqIO
from pathlib import Path
import pysam


def separate_to_files(sam_path, reads_dict, output_dir, output_format):
    output_file_handles = {}
    with pysam.AlignmentFile(sam_path, "r") as samfile:
        for record in samfile.fetch(until_eof=True):
            if record.is_unmapped or record.is_secondary or record.is_supplementary:
                continue
            if record.reference_name not in output_file_handles:
                output_path = output_dir + '/' + record.reference_name + '.' + output_format
                output_file_handle = open(output_path, 'a')
                output_file_handles[record.reference_name] = output_file_handle
            else:
                output_file_handle = output_file_handles[record.reference_name]
            SeqIO.write(reads_dict[record.query_name], output_file_handle, output_format)
    for handle in output_file_handles.values():
        handle.close()
def main():
    parser = argparse.ArgumentParser(description='Separate mapped reads into files according to mapped target.')
    parser.add_argument('-i', '--reads', type=str, required=True, help='Path to input reads.')
    parser.add_argument('-a', '--sam', type=str, required=True, help='Path to the bam/sam file.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output directory.')
    args = parser.parse_args()
    reads_path = args.reads
    sam_path = args.sam
    output_dir = args.output

    file_type = Path(reads_path).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'
    
    reads_dict = SeqIO.index(reads_path, file_type)

    separate_to_files(sam_path, reads_dict, output_dir, file_type)

if __name__ == '__main__':
    main()
