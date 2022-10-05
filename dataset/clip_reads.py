#!/usr/bin/env python
import pysam
import argparse

'''
Given a bam of reads aligned to a reference, output the reads that have the unaligned 
segments clipped off. Unmapped, seconday and supplementary alignments are skipped.
The output reads will be in the same relative strand as the alignment target.
'''
def clip_reads_fastq(sam_path):
    with pysam.AlignmentFile(sam_path, "r") as samfile:
        for record in samfile.fetch():
            if record.query_name == '927fcbe8-7891-42f8-b0d1-ac9a17b13793_2':
                print('HERE')
                print(record.is_unmapped)
            # without index will use until_eof=True, so will have unmapped
            if record.is_unmapped or record.is_secondary or record.is_supplementary:
                continue
            
            query_name = record.query_name
            aligned_segment = record.query_alignment_sequence
            quality_string = ''.join(map(lambda x: chr(x+33), record.query_alignment_qualities))

            print('@' + query_name)
            print(aligned_segment)
            print('+')
            print(quality_string)

def clip_reads_fasta(sam_path):
    with pysam.AlignmentFile(sam_path, "r") as samfile:
        for record in samfile.fetch():
            # should not have unmapped here
            if record.is_unmapped or record.is_secondary or record.is_supplementary:
                continue
            
            query_name = record.query_name
            aligned_segment = record.query_alignment_sequence
            print('>' + query_name)
            print(aligned_segment)

def main():
    parser = argparse.ArgumentParser(description='Clip the query sequences in a bam and output them as fasta/q. '
     + 'They would be in the same relative strand as the alignment target. Unmapped, seconday and supplementary are skipped.')
    parser.add_argument('-i', '--input', required=True, type=str, help='path of the input bam file.')
    parser.add_argument('-f', '--format', choices=['fasta', 'fastq'], required=True, type=str, help='Output format.')
    args = parser.parse_args()
    input_sam = args.input
    output_format = args.format

    if output_format == 'fasta':
        clip_reads_fasta(input_sam)
    elif output_format == 'fastq':
        clip_reads_fastq(input_sam)
    else:
        raise ValueError(f'Unexpected output format - {output_format}!')


if __name__ == '__main__':
    main()