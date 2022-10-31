#!/usr/bin/env python

import argparse
from Bio import SeqIO
import configparser
from dgenies_plot import startup, init_driver, plot
from haplotype_info import separate_reads_to_files
from pathlib import Path
import os
import pysam
from run_programs import run_hifiasm_homozygous_reads_only, run_pomoxis_assess_assm, run_quast, run_minimap2_reads2ref
import subprocess
import sys
import time
import matplotlib.pyplot as plt
from reads_accuracy_from_sam import estimate_error_bam_both, plot_histogram


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

def evaluate_reads_quality(reads_path, ref_path, num_threads, intermediate_dir, print_to_handle, print_to_handle_count_clip, minimap2_path='minimap2', upper_bound=100, bin_size=1):
    print('Aligning...', file=sys.stderr)
    sam_path =  intermediate_dir + '/' + Path(reads_path).stem + '_to_' + Path(ref_path).stem + '.sam'
    run_minimap2_reads2ref(reads_path, ref_path, num_threads, sam_path, minimap2_path)
        
    error_rates_no_count_clip, error_rates_count_clip = estimate_error_bam_both(sam_path, upper_bound, print_to_handle, print_to_handle_count_clip)
    
    histogram_path_count_clip = intermediate_dir + '/../' + Path(reads_path).stem + '_to_' + Path(ref_path).stem + '_count-clip.png'
    histogram_path_no_count_clip = intermediate_dir + '/../' + Path(reads_path).stem + '_to_' + Path(ref_path).stem + '_no-count-clip.png'
    plot_histogram(error_rates_count_clip, 0, upper_bound, bin_size, histogram_path_count_clip)
    plot_histogram(error_rates_no_count_clip, 0, upper_bound, bin_size, histogram_path_no_count_clip)
    

    
    
    
def assemble_hifiasm(path_to_reads, number_of_threads, output_prefix, reuse):
    return run_hifiasm_homozygous_reads_only(output_prefix, number_of_threads, [path_to_reads], reuse)

def evaluate_pomoxis(assm_paths, ref_path, output_prefixes, number_of_threads):
    for i in range(len(assm_paths)):
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
    parser = argparse.ArgumentParser(description='Evaluate the quality of haploid reads.')
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
    reads_quality_logs = dir_for_all + '/above_upper_bound.txt'
    reads_quality_logs_count_clip = dir_for_all + '/above_upper_bound_count_clip.txt'
    hifiasm_dir = dir_for_all + '/hifiasm'
    pomoxis_dir = dir_for_all + '/pomoxis'
    pomoxis_out = pomoxis_dir + '/assm'
    #pomoxis_out2 = pomoxis_dir + '/assm2'
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

    if not args.r or not os.path.isfile(reads_quality_logs):
        # Output reads quality stats
        with open(reads_quality_logs, 'w') as handle, open(reads_quality_logs_count_clip, 'w') as handle_count_clip:
            evaluate_reads_quality(
                config['DEFAULT']['reads_path'],
                config['reads_quality']['reads_ref_path'], 
                config['DEFAULT']['num_threads'],
                intermediate_dir,
                handle,
                handle_count_clip,
                config['minimap2']['path'],
                float(config['reads_quality']['upper_bound']),
                int(config['reads_quality']['bin_size'])
            )

    # Assemble with hifiasm
    # labeling of haplotype by hifiasm may not correspond to ours
    assembly_path = assemble_hifiasm(
        config['DEFAULT']['reads_path'],
        config['DEFAULT']['num_threads'],
        hifiasm_dir + '/assm',
        args.r
    )
    
    # evaluate assembly with pomoxis assess_assembly
    evaluate_pomoxis(
        [assembly_path],
        config['DEFAULT']['ref_path'],
        [pomoxis_out],
        config['DEFAULT']['num_threads']
    )

    # evaluate assembly with quast
    evaluate_quast(
        quast_dir,
        [assembly_path],
        config['DEFAULT']['ref_path'],
        config['DEFAULT']['num_threads']
    )
    
    # evaluate assembly with dgenies
    print('Starting dgenies...', file=sys.stderr)
    timeout = int(config['dgenies']['timeout']) * 3600
    try:
        with init_driver(dgenies_dir) as driver:
            dgenies_proc = startup(port_number)
            time.sleep(2)
            print('Plotting...', file=sys.stderr)
            plot(driver,
                port_number,
                config['DEFAULT']['ref_path'],
                assembly_path,
                3, timeout, 30)
            #wait_for_download(dgenies_dir, num_file1, 30)
            time.sleep(60)
    except Exception as e:
        print(e, file=sys.stderr)
    finally:
        dgenies_proc.kill()

if __name__ == '__main__':
   main()



