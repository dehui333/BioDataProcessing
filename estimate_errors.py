#!/usr/bin/env python

import argparse
import os.path
from pathlib import Path
import pysam
import subprocess
import math

'''
Dependencies:

Python libraries:
-Pysam

Tools:
-Minimap2

This script is for estimating error rates in reads w.r.t. to a reference.
(Technically an assembly can be used as 'reads' too..)

* clipped rates' percentages are w.r.t. to total length of sequences with CIGAR
  while error rates are w.r.t. to the total length of the non clipped aligned portions.

It seems that as the error rates get larger, we tend to overestimate substitution rate
and underestimate insertion and deletion rates? This has been observed in 
pomoxis scripts' assess_assembly. (could be due to misalignments at higher error rate?)


need:
1. acccept additional options
2. 

'''


'''
Given a path to a bam/sam file which contains the reads aligned to the reference, output errors rates.
Supplementary and secondary alignments are skipped.
unmapped rate calculation excludes supp and secondary.
'''
def estimate_error_bam(bam_path):
    total_aligned_len = 0 # aligned parts (not clipped)
    total_mapped_len = 0 # total len of sequences with cigar
    total_mismatch = 0
    total_ins = 0
    total_del = 0
    total_S_clip = 0
    total_H_clip = 0
    bam = pysam.AlignmentFile(bam_path) 
    num_records = 0
    num_unaligned = 0
    bins = [0 * 10]
    for alignment in bam.fetch(until_eof=True):
        if alignment.is_secondary or alignment.is_supplementary:
            continue
        num_records += 1
        if alignment.is_unmapped:
            num_unaligned += 1
            continue
        stats = alignment.get_cigar_stats()
        counts = stats[0] 
        total_mismatch += counts[8]
        total_ins += counts[1]
        total_del += counts[2]
        total_S_clip += counts[4]
        total_H_clip += counts[5]
        total_mapped_len += alignment.infer_read_length()
        total_aligned_len += alignment.query_alignment_length
        error_rate = (counts[1] + counts[2] + counts[4] + counts[5] + counts[8]) / alignment.infer_read_length()
        if (error_rate == 0):
            bins[0] += 1
            continue
        if (error_rate <= 0.05):
            bins[math.ceil(error_rate * 200) - 1] += 1
        
    bam.close()
    return total_mismatch/total_aligned_len, total_ins/total_aligned_len, total_del/total_aligned_len, total_S_clip/total_mapped_len, total_H_clip/total_mapped_len, num_unaligned/num_records, bins

'''
Given a path to a read set (fasta/q) and a reference assembly, estimate error rate in read set.
'''
def estimate_error_reads(reads_path, reference_path, keep_sam, num_threads=1): 
    reads_dir = os.path.dirname(reads_path)
    if reads_dir != '':
        reads_dir += '/' 
    sam_path =  reads_dir + Path(reads_path).stem + '_to_' + Path(reference_path).stem + '.sam'
    already_has_sam = os.path.exists(sam_path)
    if not already_has_sam:
        with open(sam_path, 'w') as sam_file:
            subprocess.run(['minimap2', '-a', '--eqx', '-t', \
                str(num_threads), reference_path, reads_path], stdout=sam_file, stderr=subprocess.DEVNULL)
        
    rates = estimate_error_bam(sam_path)
    if not keep_sam and not already_has_sam:
        subprocess.run(['rm', sam_path])
    return rates

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Estimate the error rate in a set of reads.')
    parser.add_argument('-i', '--input', type=str, help='path of the input fasta/q or bam/sam file.')
    parser.add_argument('-r', '--ref', type=str, help='path of the reference assembly.')
    parser.add_argument('-t', '--num_threads', type=int, default=1, help='number of threads for mapping. [1]')
    parser.add_argument('-k', action='store_true', help='whether to keep sam created if input is fastx.')
    parser.add_argument('-b', action='store_true', help='Input is a sam/bam.')
    args = parser.parse_args()
    if args.b == True:
        rates = estimate_error_bam(args.input)
    else:
        rates = estimate_error_reads(args.input, args.ref, args.k, args.num_threads)
    subs_rate, ins_rate, del_rate, Sc_rate, Hc_rate, unmapped_rate, bins = rates
    print(f'substitution rate: {subs_rate * 100:.3}%')
    print(f'insertion rate: {ins_rate * 100:.3}%')
    print(f'deletion rate: {del_rate * 100:.3}%')
    print(f'soft clip rate: {Sc_rate * 100:.3}%')
    print(f'hard clip rate: {Hc_rate * 100:.3}%')
    print(f'unmapped rate: {unmapped_rate * 100:.3}%')
    Sum = sum(bins)
    print(bin[0]/Sum, bin[1]/Sum, bin[2]/Sum, bin[3]/Sum, bin[4]/Sum, bin[5]/Sum, bin[6]/Sum, bin[7]/Sum, bin[8]/Sum, bin[9]/Sum)