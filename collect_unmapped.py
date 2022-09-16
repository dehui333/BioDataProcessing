#!/usr/bin/env python

import argparse
import pysam


'''
Prints unmapped reads to stdout in fasta format.
'''
def collect_unmapped(bam_path):
    bam = pysam.AlignmentFile(bam_path) 
    for record in bam.fetch(until_eof=True):
        if record.is_unmapped:
            print('>'+record.query_name)
            print(record.query_sequence)











if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collect the unmapped reads from a bam.')
    parser.add_argument('-i', '--bam', type=str, help='path of the bam/sam.')
    args = parser.parse_args()
    collect_unmapped(args.bam)