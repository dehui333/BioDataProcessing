#!/usr/bin/env python

import argparse
from collections import Counter
import matplotlib.pyplot as plt
import pysam
'''
pysam pileup seems to be slow and mem usage seems to be high. Not sure if it's my problem or what.
Process seems to get killed when the input size is large - mem use too high?

------> should use check_alignment
'''


thres = 0.2 # the second most frequently occurring char need exceed this proportion to be counted as snp
bin_size = 5000
min_coverage = 80

def find_snp_positions(bam_path, ref_path):
    ref_file = pysam.FastaFile(ref_path)
    bam_file = pysam.AlignmentFile(bam_path, "rb")
    contigs = bam_file.references
    contig_lens = bam_file.lengths
    list_of_lists_of_pos = [[] for _ in range(len(contigs))]
    contig_idx = 0
    for contig in contigs:
        if contig != '1' and contig != '2':
            continue 
        print('processing contig: ' + contig)
        for column in bam_file.pileup(contig, stepper='samtools', fastafile=ref_file, flag_filter=256+512+1024+2048):
            
            total_count = column.get_num_aligned()
            if total_count < min_coverage:
                continue
            bases = column.get_query_sequences()
            assert(total_count == len(bases))
            bases = [base.upper() for base in bases]
            counter = Counter(bases)
            if (len(counter) < 2):
                continue 
            top2 = counter.most_common(2)
            if top2[0][0] == '' or top2[1][0] == '':
                continue
            if (top2[1][1]/total_count>thres):
                list_of_lists_of_pos[contig_idx].append(column.reference_pos)
        plot_snp_positions(contigs[contig_idx], contig_lens[contig_idx], list_of_lists_of_pos[contig_idx])        
        contig_idx += 1
    ref_file.close()
    bam_file.close() 
    return contigs, contig_lens, list_of_lists_of_pos        

# should make this work on 1 contig, 1 contig len, 1 list at once
# so that some can be output faster
def plot_snp_positions(contig, contig_len, list_of_pos):
    print('plotting contig: ' + contig)
    num_bins = int(contig_len / bin_size)
    n, bins, patches = plt.hist(list_of_pos, num_bins, range=(0, contig_len))
    plt.savefig(contig + '.png')
    plt.clf()











if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count snps between two haplotypes.')
    parser.add_argument('-i', '--input', type=str, help='path of the input bam file.')
    parser.add_argument('-r', '--ref', type=str, help='path of the reference fasta file.')
    args = parser.parse_args()

    contigs, contig_lens, list_of_lists_of_pos = find_snp_positions(args.input, args.ref)
    #plot_snp_positions(contigs, contig_lens, list_of_lists_of_pos)