#!/usr/bin/env python

import argparse
import pysam

'''
Inputs: 
1. path to assembly to reference bam file with index.
2. Region, in contig:start-end format (pysam fetch format)

Outputs:
1. Regions of the assembly corresponding to the specified regions on the reference.
In format of a list of (contig, start, end) tuples 
'''
def corresponding_regions(bam_path, region):











if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Deals with regions of alignments.')
    #parser.add_argument('-i', '--assm', type=str, help='path of the assembly.')