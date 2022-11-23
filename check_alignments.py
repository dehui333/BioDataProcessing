#!/usr/bin/env python

import argparse
import edlib
from pathlib import Path
from Bio import SeqIO


class paf_record:
    def __init__(self, query_name:str, query_start:int,
                 query_end:int, same_strand:bool, target_name:str, target_start:int, target_end:int):
        self.query_name = query_name
        self.query_start = query_start
        self.query_end = query_end
        self.same_strand = same_strand
        self.target_name = target_name
        self.target_start = target_start
        self.target_end = target_end

def get_seq_file_type(seq_file_path):
    file_type = Path(seq_file_path).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'
    return file_type

def parse_paf(paf_file_path:str):
    with open(paf_file_path, 'r') as paf_file:
        for line in paf_file:
            line = line.strip().split('\t')
            
            query_name = line[0]
            query_start = int(line[2])
            query_end = int(line[3])
            strand = line[4]
            if strand == '+':
                same_strand = True
            else:
                same_strand = False 
            target_name = line[5]
            target_start = int(line[7])
            target_end = int(line[8])
        
            yield paf_record(query_name, query_start, query_end, same_strand,
                             target_name, target_start, target_end)
            
def print_alignment(target_string, query_string):
    result_dict = edlib.align(query_string, target_string, task='path')
    nice_alignment = edlib.getNiceAlignment(result_dict, query_string, target_string)
    target = nice_alignment['target_aligned']
    middle = nice_alignment['matched_aligned']
    query = nice_alignment['query_aligned']
    for i in range(0, len(target), 100):
        print(target[i:i+100])
        print(middle[i:i+100])
        print(query[i:i+100])
        print()
def main():
    parser = argparse.ArgumentParser(description='Check alignments from a paf file.')
    parser.add_argument('-i', '--reads', type=str, required=True, help='Path to input reads.')
    parser.add_argument('-r', '--read', type=str, required=True, help='A specific sequence name.')
    parser.add_argument('-p', '--paf', type=str, required=True, help='Path to paf file.')
    args = parser.parse_args()
    
    indexed_reads = SeqIO.index(args.reads, get_seq_file_type(args.reads))
    
    count = 0
    
    for record in parse_paf(args.paf):
        if record.query_name == args.read:
            count += 1
            #print(record.target_name)
            #query_string = str(indexed_reads[record.query_name].seq[record.query_start:record.query_end])
            #target_string = str(indexed_reads[record.target_name].seq[record.target_start:record.target_end])
            #print_alignment(query_string, target_string)
        if record.target_name == args.read:
            count += 1
            #print(record.query_name)
            #query_string = str(indexed_reads[record.query_name].seq[record.query_start:record.query_end])
            #target_string = str(indexed_reads[record.target_name].seq[record.target_start:record.target_end])
            #print_alignment(target_string, query_string)
    print(count)
if __name__ == '__main__':
    main()