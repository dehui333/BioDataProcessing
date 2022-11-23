#!/usr/bin/env python
import argparse
from Bio import SeqIO

def output_uniques(input_path, output_path):
    seen = set()
    with open(output_path, 'a') as output_handle:
        for record in SeqIO.parse(input_path, 'fasta'):
            if record.id not in seen:
                SeqIO.write(record, output_handle, 'fasta')
                seen.add(record.id)
            else:
                continue
                

def main():
    parser = argparse.ArgumentParser(description='Removes duplicate fasta sequences, keeping the earlier one.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to input sequences.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output path.')
    args = parser.parse_args()
    output_uniques(args.input, args.output)
    
if __name__ == '__main__':
    main()