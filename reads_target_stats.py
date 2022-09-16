#!/usr/bin/env python

import argparse
from collections import Counter
import pysam


'''
Prints the counts of where the reads are mapped to.
'''
def count_targets(bam_path):
    bam = pysam.AlignmentFile(bam_path) 
    num_unmapped = 0
    primary_counter = Counter()
    secondary_counter = Counter()
    supplementary_counter = Counter()

    def update_target_count(target, counter):
        counter.update({target: 1})

    for record in bam.fetch(until_eof=True):
        if record.is_unmapped:
            num_unmapped +=1
            continue
        if record.is_secondary:
            update_target_count(record.reference_name, secondary_counter)
            continue
        if record.is_supplementary:
            update_target_count(record.reference_name, supplementary_counter)
            continue
        update_target_count(record.reference_name, primary_counter)
    print('num unmapped: ', num_unmapped)
    print('primary alignments: ')
    print_counts(primary_counter)
    print('secondary alignments: ')
    print_counts(secondary_counter)
    print('supplementary alignments: ')
    print_counts(supplementary_counter)

def print_counts(counter):
    for target, count in counter.items():
        print(target, ': ', count)










if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count the occurrences of where reads are mapped to.')
    parser.add_argument('-i', '--bam', type=str, help='path of the bam/sam.')
    args = parser.parse_args()
    count_targets(args.bam)