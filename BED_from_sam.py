#!/usr/bin/env python

import argparse
import pysam
import sys

'''
Takes in an alignment file (sam/bam).
Extract the regions with (primary) alignments into a BED file.

BED output is printed to stdout
number of unmapped queries is printed to stderr

'''

def print_aligned_regions(sam_path, out_handle, err_handle):
    num_unmapped = 0
    with pysam.AlignmentFile(sam_path) as sam:
        for record in sam.fetch(until_eof=True):
            if record.is_secondary or record.is_supplementary:
                continue
            if record.is_unmapped:
                num_unmapped += 1
                continue
            print(record.reference_name + '\t' + str(record.reference_start) + '\t' + str(record.reference_end), file=out_handle)
    print('num unmapped ', num_unmapped, file=err_handle)

def main():
    parser = argparse.ArgumentParser(description='Output the aligned regions on the reference in an alignment file in BED format.'
                                    + ' BED output is printed to stdout.'
                                    + ' Number of unmapped queries is printed to stderr.')
    parser.add_argument('-i', '--input', type=str, help='path of the input bam/sam file.')
    args = parser.parse_args()
    print_aligned_regions(args.input, sys.stdout, sys.stderr)

if __name__ == '__main__':
    main()