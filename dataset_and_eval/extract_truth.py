#!/usr/bin/env python

import argparse
import pysam


'''
Extract the aligned segments on the reference as truth reads.
'''

def extract(reads2ref_bam_path):  
    with pysam.AlignmentFile(reads2ref_bam_path) as bam:
        for record in bam.fetch():
            if record.is_unmapped or record.is_secondary or record.is_supplementary:
                continue
            print('>'+record.query_name)
            print(record.get_reference_sequence().upper())

def main():
    parser = argparse.ArgumentParser(description='Extract the truth reads from reads2referece alignment.')
    parser.add_argument('-i', '--input', type=str, help='Path to the bam/sam.')
    args = parser.parse_args()

    extract(args.input)

if __name__ == '__main__':
    main()