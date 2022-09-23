#!/usr/bin/env python

import argparse
import os.path
from pathlib import Path
import pysam


'''
Extract the aligned segments on the reference as truth reads.
'''

def extract(reads2ref_bam_path):
    bam = pysam.AlignmentFile(reads2ref_bam_path) 
    for record in bam.fetch(until_eof=True):
        if record.is_unmapped:
            continue
        if record.is_secondary:
            continue
        if record.is_supplementary:
            continue
        
        print('>'+record.query_name + ' ')
        print(record.get_reference_sequence())




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract the truth reads from reads2trimmed alignment.')
    parser.add_argument('-i', '--input', type=str, help='Path to the bam.')
    args = parser.parse_args()

    extract(args.input)