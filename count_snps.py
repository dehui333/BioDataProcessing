#!/usr/bin/env python

import argparse
from Bio import SeqIO
from check_alignment import get_snp_pos
import matplotlib.pyplot as plt


'''
pysam pileup seems to be slow and mem usage seems to be high. Not sure if it's my problem or what.
Process seems to get killed when the input size is large - mem use too high?

------> should use check_alignment
'''


thres = 0.3 # the second most frequently occurring char need exceed this proportion to be counted as snp
bin_size = 5000
min_coverage = 80
min_mapq = 10
skip_HP = 1
skip_indel = 0

def plot_snp_positions(bam_path, ref_path):
    contigs_dict = SeqIO.index(ref_path, 'fasta')
    contig_names = []
    contig_lens = []
    list_of_lists_of_pos = []
    for contig_name, contig in contigs_dict.items():
        contig = str(contig.seq)
        contig_len = len(contig)
        positions = get_snp_pos(args.input, contig_name, contig, contig_len, min_mapq, min_coverage, thres, skip_HP, skip_indel)
        plot(contig_name, contig_len, positions)
        contig_names.append(contig_name)
        contig_lens.append(contig_len)
        list_of_lists_of_pos.append(positions)
    return contig_names, contig_lens, list_of_lists_of_pos        

# should make this work on 1 contig, 1 contig len, 1 list at once
# so that some can be output faster
def plot(contig_name, contig_len, list_of_pos):
    print('plotting contig: ' + contig_name)
    num_bins = int(contig_len / bin_size)
    n, bins, patches = plt.hist(list_of_pos, num_bins, range=(0, contig_len))
    plt.savefig(args.output + '/' + contig_name + '.png')
    plt.clf()







'''
 if (!PyArg_ParseTuple(args, "sssllldi", &bam_path, &contig_name, &contig,
                          &contig_len, &min_mapq, &min_coverage, &min_alt_prop, &skip_HP))
'''



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count snps between two haplotypes.')
    parser.add_argument('-i', '--input', type=str, help='path of the input bam file.')
    parser.add_argument('-r', '--ref', type=str, help='path of the reference fasta file.')
    parser.add_argument('-o', '--output', type=str, help='path of the output directory.')
    args = parser.parse_args()


    plot_snp_positions(args.input, args.ref)
    