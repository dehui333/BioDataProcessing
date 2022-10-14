#!/usr/bin/env python

import argparse
from Bio import SeqIO
import configparser
from dgenies_plot import startup, wait_for_download, init_driver, plot
from haplotype_info import separate_reads_to_files
from pathlib import Path
import os
import pysam
from run_programs import run_hifiasm_hetero_reads_only, run_pomoxis_assess_assm, run_quast, run_minimap2_reads2ref
import subprocess
import sys
import time

'''
Other dependencies(put in PATH):
1. minimap2 (can specify in config file for evaluate.py)
2. quast
3. hifiasm


The config file for dgenies and the its minimap2 instance is at the {path of dgenies.__file__}/../etc/genies
(to change thread number, max size of upload file etc)
The directories need to be absolute path

***** Should try to handle exceptions etc more cleanly..
***** The evaluation on assemblies can be refactored and put together as an independent functionality
'''

  
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

def evaluate_reads_quality(reads_path, ref1_path, ref2_path, name1, name2, num_threads, intermediate_dir, print_to_handle=sys.stdout, minimap2_path='minimap2'):
    print('Separating reads...', file=sys.stderr)
    reads1, reads2 = separate_reads_to_files(reads_path, name1, name2, intermediate_dir)
    print('Aligning...', file=sys.stderr)
    sam1_path =  intermediate_dir + '/' + Path(reads1).stem + '_to_' + Path(ref1_path).stem + '.sam'
    sam2_path =  intermediate_dir + '/' + Path(reads2).stem + '_to_' + Path(ref2_path).stem + '.sam'
    run_minimap2_reads2ref(reads1, ref1_path, num_threads, sam1_path, minimap2_path)
    run_minimap2_reads2ref(reads2, ref2_path, num_threads, sam2_path, minimap2_path)
    stats1 = basic_acc_stats(sam1_path)
    stats2 = basic_acc_stats(sam2_path)
    print(f'For {name1}:', file=print_to_handle)
    print_basic_acc_stats(stats1, print_to_handle=print_to_handle)
    print(f'For {name1}:', file=print_to_handle)
    print_basic_acc_stats(stats2, print_to_handle=print_to_handle)

def assemble_hifiasm(path_to_reads, number_of_threads, output_prefix, reuse):
    return run_hifiasm_hetero_reads_only(output_prefix, number_of_threads, [path_to_reads], reuse)

def evaluate_pomoxis(assm_paths, ref_path, output_prefixes, number_of_threads):
    for i in range(2):
        run_pomoxis_assess_assm(assm_paths[i], ref_path, number_of_threads, output_prefixes[i])
    
def evaluate_quast(output_dir, assm_paths, ref_path, num_threads):
    run_quast(output_dir, assm_paths, ref_path, num_threads)

# from https://stackoverflow.com/questions/2470971/fast-way-to-test-if-a-port-is-in-use-using-python
def is_port_in_use(port: int) -> bool:
    import socket
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(('localhost', port)) == 0

def main():
     # Parsing
    parser = argparse.ArgumentParser(description='Evaluate the quality of diploid reads.')
    parser.add_argument('-c', '--config', type=str, help='Path to config file.')
    parser.add_argument('-r', action='store_true', help='Try to reuse existing files.')
    
    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(args.config)

    # checking port for dgenies
    port_number = int(config['dgenies']['port_number'])
    if is_port_in_use(port_number):
        print(f'Port {port_number} is in use!', file=sys.stderr)
        exit(1) 
    
    # Prepare output directories
    dir_for_all = config['DEFAULT']['output_dir']


    # outputs
    reads_quality_stats = dir_for_all + '/reads_quality.txt'
    hifiasm_dir = dir_for_all + '/hifiasm'
    pomoxis_dir = dir_for_all + '/pomoxis'
    pomoxis_out1 = pomoxis_dir + '/assm1'
    pomoxis_out2 = pomoxis_dir + '/assm2'
    intermediate_dir = dir_for_all + '/intermediates'
    dgenies_dir = dir_for_all + '/dgenies'
    dgenies_log = dgenies_dir + '/log.txt'
    quast_dir = dir_for_all + '/quast' # will be created by quast
    
    if not dir_for_all:
        print(f'Output directory not specified!', file=sys.stderr)
        exit(1)
    elif os.path.isdir(dir_for_all):
        if not args.r:
            print(f'{dir_for_all} already exists!', file=sys.stderr)
            exit(1)
    else:
        os.mkdir(dir_for_all)
        os.mkdir(intermediate_dir)
        os.mkdir(pomoxis_dir)
        os.mkdir(hifiasm_dir)
        os.mkdir(dgenies_dir)

    if not args.r or not os.path.isfile(reads_quality_stats):
        # Output reads quality stats
        with open(reads_quality_stats, 'w') as handle:
            evaluate_reads_quality(
                config['DEFAULT']['reads_path'],
                config['reads_quality']['ref1_path'],
                config['reads_quality']['ref2_path'],
                config['reads_quality']['name1'],
                config['reads_quality']['name2'], 
                config['DEFAULT']['num_threads'],
                intermediate_dir,
                handle,
                config['minimap2']['path']
            )

    # Assemble with hifiasm
    # labeling of haplotype by hifiasm may not correspond to ours
    assembly_paths = list(assemble_hifiasm(
        config['DEFAULT']['reads_path'],
        config['DEFAULT']['num_threads'],
        hifiasm_dir + '/assm',
        args.r
    ))
    
    # evaluate assembly with pomoxis assess_assembly
    evaluate_pomoxis(
        assembly_paths,
        config['DEFAULT']['ref_path'],
        [pomoxis_out1, pomoxis_out2],
        config['DEFAULT']['num_threads']
    )

    # evaluate assembly with quast
    evaluate_quast(
        quast_dir,
        assembly_paths,
        config['DEFAULT']['ref_path'],
        config['DEFAULT']['num_threads']
    )
    
    # evaluate assembly with dgenies
    print('Starting dgenies...', file=sys.stderr)
    timeout = config['dgenies']['timeout'] * 3600
    try:
        with init_driver(dgenies_dir) as driver:
            dgenies_proc = startup(port_number)
            time.sleep(2)
            print('Plotting...', file=sys.stderr)
            plot(driver,
                port_number,
                config['DEFAULT']['ref_path'],
                assembly_paths[0],
                3, timeout, 30)
            #wait_for_download(dgenies_dir, num_file1, 30)
            time.sleep(20)
            print('Plotting...', file=sys.stderr)
            plot(driver,
                port_number,
                config['DEFAULT']['ref_path'],
                assembly_paths[1],
                3, timeout, 30)
            #wait_for_download(dgenies_dir, num_file1+num_file2, 30)
            time.sleep(60)
    except Exception as e:
        print(e, file=sys.stderr)
    finally:
        dgenies_proc.kill()

if __name__ == '__main__':
   main()





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


Separate reads into different files based on haplotype, return the file names.

def separate_reads_to_files(reads_path, key, tag1, tag2, output_dir):
    file_type = Path(reads_path).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'

    out1 = tag1 + '.' + file_type
    out2 = tag2 + '.' + file_type
    out1 = output_dir + '/' + Path(reads_path).stem + '-' + out1
    out2 = output_dir + '/' + Path(reads_path).stem + '-' + out2
    if os.path.isfile(out1) and os.path.isfile(out2):
        print(f'{out1} and {out2} already exists. Using existing ones.', file=sys.stderr)
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

Align reads to reference then return the path to the sam.
maybe should move to run_programs.py

def align_reads2ref(reads_path, ref_path, num_threads, output_dir):
    #reads_dir = os.path.dirname(reads_path)
    #if reads_dir != '':
    #    reads_dir += '/' 
    sam_path =  output_dir + Path(reads_path).stem + '_to_' + Path(ref_path).stem + '.sam'
    sam_path = sam_path
    already_has_sam = os.path.exists(sam_path)
    if not already_has_sam:
        with open(sam_path, 'w') as sam_file:
            subprocess.run(['minimap2', '-a', '--eqx', '-t', \
                str(num_threads), ref_path, reads_path], stdout=sam_file)
    else:
        print(sam_path + ' already exists! Using existing one.', file=sys.stderr)
    return sam_path
'''      