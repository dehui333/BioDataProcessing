#!/usr/bin/env python

import argparse
from unittest import result
from Bio import SeqIO
import edlib
from pathlib import Path
import matplotlib.pyplot as plt
import sys


def get_error_rates(seq_file_path, truth_dict, upper_bound, logfile_handle, print_alignment=False):
    num_above = 0
    num_unmapped = 0
    error_rates = []
    task = 'distance'
    if print_alignment:
        task = 'path'
    for record in SeqIO.parse(seq_file_path, get_seq_file_type(seq_file_path)):
        if record.id not in truth_dict:
            num_unmapped += 1
            continue
        query_string = str(record.seq)
        truth_string = str(truth_dict[record.id].seq)
        result_dict = edlib.align(str(record.seq), truth_string, task=task)
        error_rate = result_dict['editDistance']/len(query_string) * 100
     
        if error_rate < upper_bound:
            error_rates.append(error_rate)
        else:
            num_above+=1
            print(record.id, error_rate, file=logfile_handle)
            if print_alignment:
                nice_alignment = edlib.getNiceAlignment(result_dict, query_string, truth_string)
                target = nice_alignment['target_aligned']
                middle = nice_alignment['matched_aligned']
                query = nice_alignment['query_aligned']
                for i in range(0, len(target), 100):
                    print(target[i:i+100], file=logfile_handle)
                    print(middle[i:i+100], file=logfile_handle)
                    print(query[i:i+100], file=logfile_handle)
                    print(file=logfile_handle)
    print('num above upper bound' , num_above, file=logfile_handle)
    print('num unmapped', num_unmapped, file=logfile_handle)
    return error_rates

def get_seq_file_type(seq_file_path):
    file_type = Path(seq_file_path).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'
    return file_type

def plot_histogram(values, range_start, range_end, bin_size, output_path):
    num_bins = int((range_end- range_start) / bin_size)
    n, bins, patches = plt.hist(values, num_bins, range=(range_start, range_end))
    plt.savefig(output_path)
    plt.clf()

def main():
    parser = argparse.ArgumentParser(description='Evaluate reads accuracy w.r.t. a set of ground-truth.')
    parser.add_argument('-i', '--reads', type=str, required=True, help='Path to input reads.')
    parser.add_argument('-r', '--truth', type=str, required=True, help='Path to ground-truth reads.')
    parser.add_argument('-u', '--upper', type=float, required=True, help='The upper bound percentage of error. Those above will not be plotted in histogram.')
    parser.add_argument('-b', '--bin', type=float, required=True, help='The bin size in percentages.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output path for histogram.')
    parser.add_argument('-p', '--print', action='store_true', help='Print alignments that have higher than upper bound error rate to stdout.')
    args = parser.parse_args()
    reads_path = args.reads
    truth_path = args.truth
    out_path = args.output
    error_upper_bound = args.upper
    truth_dict = SeqIO.index(truth_path, get_seq_file_type(truth_path))
    # in percentages
    error_rates = get_error_rates(reads_path, truth_dict, error_upper_bound, sys.stdout, args.print)
    plot_histogram(error_rates, 0, error_upper_bound, args.bin, out_path)
    
if __name__ == '__main__':
    main()
    
    

    
