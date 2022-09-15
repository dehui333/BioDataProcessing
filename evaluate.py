#!/usr/bin/env python

import argparse
from Bio import SeqIO
from pathlib import Path
import os.path
import pysam
import subprocess

PREFIX = 'EVAL_'

'''
To evaluate the quality of reads w.r.t. alignment to trimmed assembly.

1. put reads from different strains into different files using info in the description.
2. align each read set onto the trimmed assembly of that strain.
3. calculate accuracy.
'''
# Return the haplotype number (1 or 2) by looking at the description of a sequence.
def get_hap(description):
    key_start = description.find(args.kv[0]+'=')
    if key_start == -1:
        print('Sequence has no "' + args.kv[0] + '" field in description!')
        exit(0)
    value_start = key_start + len(args.kv[0]) + 1
    if description.find(args.kv[1], value_start) != -1:
        return 1
    elif description.find(args.kv[2], value_start) != -1:
        return 2
    print('Sequence does not follow haplotypes as specified.')
    exit(0)

'''
Separate reads into different files based on haplotype, return the file names.
'''
def separate_reads_to_files(reads_path):
    file_type = Path(reads_path).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'

    out1 = args.kv[1] + '.' + file_type
    out2 = args.kv[2] + '.' + file_type
    out1 = PREFIX + out1
    out2 = PREFIX + out2
    if os.path.isfile(out1) and os.path.isfile(out2):
        print(f'{out1} or {out2} already exists. Using existing one.')
        return out1, out2
    elif os.path.isfile(out1) or os.path.isfile(out2):
        print(f'One of {out1} and {out2} already exists!')
        exit(0)
    with open(out1, 'a') as h1_out, open(out2, 'a') as h2_out:
        for record in SeqIO.parse(reads_path, file_type):
            hap_number = get_hap(record.description)
            if hap_number == 1:
                SeqIO.write(record, h1_out, file_type)
            elif hap_number == 2:
                SeqIO.write(record, h2_out, file_type)
            else:
                print('??')
                exit(0)
    return out1, out2
'''
Align reads to reference then return the path to the sam.
'''
def align_reads2ref(reads_path, ref_path):
    #reads_dir = os.path.dirname(reads_path)
    #if reads_dir != '':
    #    reads_dir += '/' 
    sam_path =  Path(reads_path).stem + '_to_' + Path(ref_path).stem + '.sam'
    sam_path = PREFIX + sam_path
    already_has_sam = os.path.exists(sam_path)
    if not already_has_sam:
        with open(sam_path, 'w') as sam_file:
            subprocess.run(['minimap2', '-a', '--eqx', '-t', \
                str(args.threads), ref_path, reads_path], stdout=sam_file)
    else:
        print(sam_path + ' already exists! Using existing one.')
    return sam_path
        
def basic_acc_stats(bam_path):
    total_aligned_len = 0 # aligned parts (not clipped)
    total_mapped_len = 0 # total len of sequences with cigar
    total_mismatch = 0
    total_ins = 0
    total_del = 0
    total_S_clip = 0
    total_H_clip = 0
    bam = pysam.AlignmentFile(bam_path) 
    num_records = 0 # exclude secondary/supplementary
    num_unaligned = 0
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
        
    bam.close()
    return total_mismatch, total_ins, total_del, total_S_clip, total_H_clip, total_aligned_len, total_mapped_len, num_unaligned, num_records

def print_basic_acc_stats(stats):
    total_mismatch = stats[0]
    total_ins = stats[1]
    total_del = stats[2]
    total_S_clip = stats[3]
    total_H_clip = stats[4]
    total_aligned_len = stats[5]
    total_mapped_len = stats[6]
    num_unaligned = stats[7]
    num_records = stats[8]
    print(f'Percentage mismatch: {total_mismatch/total_aligned_len * 100:.3}%')
    print(f'Percentage ins: {total_ins/total_aligned_len * 100:.3}%')
    print(f'Percentage del: {total_del/total_aligned_len * 100:.3}%')
    print(f'Percentage softclipped: {total_S_clip/total_mapped_len * 100:.3}%')
    print(f'Percentage hardclipped: {total_H_clip/total_mapped_len * 100:.3}%')
    print(f'Percentage unmapped: {num_unaligned/num_records * 100:.3}%')

def evaluate_reads_quality():
    print('Separating reads...')
    reads1, reads2 = separate_reads_to_files(args.reads)
    print('Aligning...')
    sam1 = align_reads2ref(reads1, args.ref1)
    sam2 = align_reads2ref(reads2, args.ref2)
    stats1 = basic_acc_stats(sam1)
    stats2 = basic_acc_stats(sam2)
    print(f'For {args.kv[1]}:')
    print_basic_acc_stats(stats1)
    print(f'For {args.kv[2]}:')
    print_basic_acc_stats(stats2)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Evaluate the quality of reads w.r.t. alignment to trimmed assembly.')
    parser.add_argument('-i', '--reads', type=str, help='Path to the reads.', required=True)
    parser.add_argument('--kv', nargs=3, help='Determines which key determines haplotype and which value map to which haplotype. e.g. strain S288C FSY1742. The two values should not be prefix of each other.', required=True)
    parser.add_argument('--ref1', type=str, help='Path to the reference (trimmed assembly) of haplotype1.', required=True)
    parser.add_argument('--ref2', type=str, help='Path to the reference (trimmed assembly) of haplotype2.', required=True)
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use.')
    args = parser.parse_args()
    evaluate_reads_quality()