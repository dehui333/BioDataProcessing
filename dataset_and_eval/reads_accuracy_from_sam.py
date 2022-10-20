#!/usr/bin/env python

import argparse
import os.path
from pathlib import Path
import pysam
import matplotlib.pyplot as plt
import sys

# should make version to get both errors including/excluding clipping in one function
def estimate_error_bam(bam_path, count_clip, upper_bound, logfile_handle):
    error_rates = []
    num_unmapped = 0
    num_above_bound = 0
    with pysam.AlignmentFile(bam_path) as bam:
        for record in bam.fetch(until_eof=True): # if not index, same as until_eof=True
            if record.is_secondary or record.is_supplementary:
                continue
            if record.is_unmapped:
                num_unmapped += 1
                continue
            stats = record.get_cigar_stats()
            counts = stats[0] 
            num_mismatch = counts[8]
            num_ins = counts[1]
            num_del = counts[2]
            num_soft_clip = counts[4]
            num_hard_clip = counts[5]
            total_query_len = record.infer_read_length()
            aligned_query_len = record.query_alignment_length
            if count_clip:
                error_rate = (num_mismatch + num_ins + num_del \
                + num_soft_clip + num_hard_clip) / total_query_len * 100
            else:
                error_rate = (num_mismatch + num_ins + num_del) / aligned_query_len * 100
            if error_rate < upper_bound:
                error_rates.append(error_rate)
            else:
                print(record.query_name, error_rate, file=logfile_handle) 
                num_above_bound += 1
    print('num above bound ', num_above_bound, file=logfile_handle)
    print('num unmapped ', num_unmapped, file=logfile_handle) 
    return error_rates

# version to get both errors including/excluding clipping in one function
def estimate_error_bam_both(bam_path, upper_bound, logfile_handle, logfile_handle_count_clip):
    error_rates = []
    error_rates_count_clip = []
    num_unmapped = 0
    num_above_bound = 0
    num_above_bound_count_clip = 0
    with pysam.AlignmentFile(bam_path) as bam:
        for record in bam.fetch(until_eof=True): # if not index, same as until_eof=True
            if record.is_secondary or record.is_supplementary:
                continue
            if record.is_unmapped:
                num_unmapped += 1
                continue
            stats = record.get_cigar_stats()
            counts = stats[0] 
            num_mismatch = counts[8]
            num_ins = counts[1]
            num_del = counts[2]
            num_soft_clip = counts[4]
            num_hard_clip = counts[5]
            total_query_len = record.infer_read_length()
            aligned_query_len = record.query_alignment_length
            
            error_rate_count_clip = (num_mismatch + num_ins + num_del \
            + num_soft_clip + num_hard_clip) / total_query_len * 100
    
            error_rate = (num_mismatch + num_ins + num_del) / aligned_query_len * 100
            if error_rate < upper_bound:
                error_rates.append(error_rate)
            else:
                print(record.query_name, error_rate, file=logfile_handle) 
                num_above_bound += 1
            
            if error_rate_count_clip < upper_bound:
                error_rates_count_clip.append(error_rate_count_clip)
            else:
                print(record.query_name, error_rate_count_clip, file=logfile_handle_count_clip)
                num_above_bound_count_clip += 1
    print('num above bound ', num_above_bound, file=logfile_handle)
    print('num unmapped ', num_unmapped, file=logfile_handle) 
    
    print('num above bound ', num_above_bound_count_clip, file=logfile_handle_count_clip)
    print('num unmapped ', num_unmapped, file=logfile_handle_count_clip)
    return error_rates, error_rates_count_clip

def plot_histogram(values, range_start, range_end, bin_size, output_path):
    num_bins = int((range_end- range_start) / bin_size)
    n, bins, patches = plt.hist(values, num_bins, range=(range_start, range_end))
    plt.savefig(output_path)
    plt.clf()

def main():
    parser = argparse.ArgumentParser(description='Estimate the error rate in a set of reads from its alignment to a reference.')
    parser.add_argument('-i', '--input', type=str, help='path of the input bam/sam file.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output path for histogram.')
    parser.add_argument('-b', '--bin', type=float, required=True, help='Bin size, in percentage.')
    parser.add_argument('-c', '--clip', action='store_true', help='Count clipping as errors.')
    parser.add_argument('-u', '--upper', type=float, required=True, help='The upper bound percentage of error. Those above will not be plotted in histogram and their info is logged.')
    args = parser.parse_args()
    error_rates = estimate_error_bam(args.input, args.clip, args.upper, sys.stdout)
    plot_histogram(error_rates, 0, args.upper, args.bin, args.output)

if __name__ == '__main__':
    main()
    