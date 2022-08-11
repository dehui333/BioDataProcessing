from add_errors import add_errors
import argparse
import configparser
from simulate_perfect_reads import seqreq_perfect_reads

'''
This script generates simulated reads with simple errors from an input genome fasta file.
Uses the same config file format as simulating perfect reads, though its output path will be ignored.
'''



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads with simple errors from a genome fasta file.')
    # for simulating reads
    parser.add_argument('-c', '--config', type=str, help='path to the config file for settings on simulating reads.')
    
    # for adding error
    parser.add_argument('-p', '--perfect_output', type=str, help='output path of the perfect reads.')
    parser.add_argument('-e', '--error_output', type=str, help='output path of the reads with error.')
    parser.add_argument('--sub', type=float, default=0.03, help='substitution rate.')
    parser.add_argument('--ins', type=float, default=0.03, help='insertion rate.')
    parser.add_argument('--dele', type=float, default=0.03, help='deletion rate.')
    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(args.config)
    
    # force perfect reads output path to follow arguments 
    config['seqrequester']['output_path'] = args.perfect_output
    

    seqreq_perfect_reads(config)
    add_errors(args.perfect_output, args.error_output, args.sub, args.ins, args.dele)

