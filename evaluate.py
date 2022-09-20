#!/usr/bin/env python

import argparse
from Bio import SeqIO
import configparser
from pathlib import Path
import os
import pysam
from run_programs import run_hifiasm_hetero_reads_only, run_pomoxis_assess_assm, run_quast
import subprocess
import sys


'''
To evaluate the quality of reads w.r.t. alignment to trimmed assembly.

1. put reads from different strains into different files using info in the description.
2. align each read set onto the trimmed assembly of that strain.
3. calculate accuracy.
'''
# Return the haplotype number (1 or 2) by looking at the description of a sequence.
def get_hap(description, key, tag1, tag2):
    key_start = description.find(key+'=')
    if key_start == -1:
        print('Sequence has no "' + key + '" field in description!', file=sys.stderr)
        exit(1)
    value_start = key_start + len(key) + 1
    if description.find(tag1, value_start) != -1:
        return 1
    elif description.find(tag2, value_start) != -1:
        return 2
    print('Sequence does not follow haplotypes as specified.', file=sys.stderr)
    exit(1)

'''
Separate reads into different files based on haplotype, return the file names.
'''
def separate_reads_to_files(reads_path, key, tag1, tag2):
    file_type = Path(reads_path).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'

    out1 = tag1 + '.' + file_type
    out2 = tag2 + '.' + file_type
    out1 = INTERMEDIATE_PREFIX + Path(reads_path).stem + '-' + out1
    out2 = INTERMEDIATE_PREFIX + Path(reads_path).stem + '-' + out2
    if os.path.isfile(out1) and os.path.isfile(out2):
        print(f'{out1} or {out2} already exists. Using existing one.', file=sys.stderr)
        return out1, out2
    elif os.path.isfile(out1) or os.path.isfile(out2):
        print(f'One of {out1} and {out2} already exists!', file=sys.stderr)
        exit(1)
    with open(out1, 'a') as h1_out, open(out2, 'a') as h2_out:
        for record in SeqIO.parse(reads_path, file_type):
            hap_number = get_hap(record.description, key, tag1, tag2)
            if hap_number == 1:
                SeqIO.write(record, h1_out, file_type)
            elif hap_number == 2:
                SeqIO.write(record, h2_out, file_type)
            else:
                print('??', file=sys.stderr)
                exit(1)
    return out1, out2
'''
Align reads to reference then return the path to the sam.
'''
def align_reads2ref(reads_path, ref_path, num_threads):
    #reads_dir = os.path.dirname(reads_path)
    #if reads_dir != '':
    #    reads_dir += '/' 
    sam_path =  INTERMEDIATE_PREFIX + Path(reads_path).stem + '_to_' + Path(ref_path).stem + '.sam'
    sam_path = sam_path
    already_has_sam = os.path.exists(sam_path)
    if not already_has_sam:
        with open(sam_path, 'w') as sam_file:
            subprocess.run(['minimap2', '-a', '--eqx', '-t', \
                str(num_threads), ref_path, reads_path], stdout=sam_file)
    else:
        print(sam_path + ' already exists! Using existing one.', file=sys.stderr)
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

def print_basic_acc_stats(stats, print_to_handle):
    total_mismatch = stats[0]
    total_ins = stats[1]
    total_del = stats[2]
    total_S_clip = stats[3]
    total_H_clip = stats[4]
    total_aligned_len = stats[5]
    total_mapped_len = stats[6]
    num_unaligned = stats[7]
    num_records = stats[8]
    print(f'Percentage mismatch: {total_mismatch/total_aligned_len * 100:.3}%', file=print_to_handle)
    print(f'Percentage ins: {total_ins/total_aligned_len * 100:.3}%', file=print_to_handle)
    print(f'Percentage del: {total_del/total_aligned_len * 100:.3}%', file=print_to_handle)
    print(f'Percentage softclipped: {total_S_clip/total_mapped_len * 100:.3}%', file=print_to_handle)
    print(f'Percentage hardclipped: {total_H_clip/total_mapped_len * 100:.3}%', file=print_to_handle)
    print(f'Percentage unmapped: {num_unaligned/num_records * 100:.3}%', file=print_to_handle)

def evaluate_reads_quality(reads_path, ref1_path, ref2_path, key, tag1, tag2, num_threads, print_to_handle=sys.stdout):
    print('Separating reads...', file=sys.stderr)
    reads1, reads2 = separate_reads_to_files(reads_path, key, tag1, tag2)
    print('Aligning...', file=sys.stderr)
    sam1 = align_reads2ref(reads1, ref1_path, num_threads)
    sam2 = align_reads2ref(reads2, ref2_path, num_threads)
    stats1 = basic_acc_stats(sam1)
    stats2 = basic_acc_stats(sam2)
    print(f'For {tag1}:', file=print_to_handle)
    print_basic_acc_stats(stats1, print_to_handle=print_to_handle)
    print(f'For {tag2}:', file=print_to_handle)
    print_basic_acc_stats(stats2, print_to_handle=print_to_handle)

def assemble_hifiasm(path_to_reads, number_of_threads, output_prefix):
    print('Assembling with hifiasm...', file=sys.stderr)
    return run_hifiasm_hetero_reads_only(output_prefix, number_of_threads, [path_to_reads])

def evaluate_pomoxis(assm_paths, ref_path, output_prefixes, number_of_threads):
    for i in range(2):
        run_pomoxis_assess_assm(assm_paths[i], ref_path, number_of_threads, output_prefixes[i])
    
def evaluate_quast(output_dir, assm_paths, ref_path, num_threads):
    run_quast(output_dir, assm_paths, ref_path, num_threads)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Evaluate the quality of diploid reads.')
    parser.add_argument('-c', '--config', type=str, help='Path to config file.')
    
    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(args.config)
    
    dir_for_all = config['DEFAULT']['output_dir']
    
    if not dir_for_all:
        print(f'Output directory not specified!', file=sys.stderr)
        exit(1)
    elif os.path.isdir(dir_for_all):
        print(f'{dir_for_all} already exists!', file=sys.stderr)
        exit(1)
    else:
        os.mkdir(dir_for_all)
        global INTERMEDIATE_PREFIX
        INTERMEDIATE_PREFIX = dir_for_all + '/intermediates/' 
        os.mkdir(INTERMEDIATE_PREFIX)
        os.mkdir(dir_for_all + '/pomoxis')
        os.mkdir(dir_for_all + '/hifiasm')
    

    with open(dir_for_all + '/reads_quality.txt', 'w') as handle:
        evaluate_reads_quality(
            config['DEFAULT']['reads_path'],
            config['reads_quality']['ref1_path'],
            config['reads_quality']['ref2_path'],
            config['reads_quality']['key'],
            config['reads_quality']['tag1'],
            config['reads_quality']['tag2'], 
            config['DEFAULT']['num_threads'],
            handle
        )

    # labeling of haplotype by hifiasm may not correspond to ours
    assembly_paths = list(assemble_hifiasm(
        config['DEFAULT']['reads_path'],
        config['DEFAULT']['num_threads'],
        dir_for_all + '/hifiasm/assm'
    ))
    
    
    evaluate_pomoxis(
        assembly_paths,
        config['DEFAULT']['ref_path'],
        [dir_for_all + '/pomoxis/assm1', dir_for_all + '/pomoxis/assm2'],
        config['DEFAULT']['num_threads']
    )


    evaluate_quast(
        dir_for_all + '/quast',
        assembly_paths,
        config['DEFAULT']['ref_path'],
        config['DEFAULT']['num_threads']
    )

    