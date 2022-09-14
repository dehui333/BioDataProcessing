#!/usr/bin/env python

import argparse


'''
To evaluate the quality of reads w.r.t. alignment to trimmed assembly.

1. put reads from different strains into different files using info in the description.
2. align each read set onto the trimmed assembly of that strain.
3. calculate accuracy.
'''


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Evaluate the quality of reads w.r.t. alignment to trimmed assembly.')
    parser.add_argument('-i', '--reads', type=str, help='Path to the reads.', required=True)
    parser.add_argument('--kv', nargs=3, help='Determines which key determines haplotype and which value map to which haplotype. e.g. strain S288C FSY1742', required=True)
    parser.add_argument('--ref1', type=str, help='Path to the reference (trimmed assembly) of haplotype1.', required=True)
    parser.add_argument('--ref2', type=str, help='Path to the reference (trimmed assembly) of haplotype2.', required=True)
    args = parser.parse_args()
    print(args.kv)
    #evaluate_quality(args.reads, args.bam)