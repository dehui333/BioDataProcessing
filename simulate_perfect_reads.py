import argparse
from Bio import SeqIO
from collections import Counter
import configparser
import os.path
import subprocess

'''
Dependencies:

Python libraries:
-BioPython

Tools:
-Seqrequester

This is a script for generating perfect simulated reads using existing tools, currently: seqrequester.

'''
def get_genome_length(genome_path):
    total_len = 0
    for record in SeqIO.parse(genome_path, 'fasta'):
        total_len += len(record)
    return total_len

# For Seqrequester
def get_length_profile(reads_path):
    profile_path = reads_path + '.profile'
    if os.path.exists(profile_path):
        print(f'Using {profile_path} as reads length distribution profile...')
        return profile_path
    print('Generating reads length distribution profile...')
    lengths = []
    for record in SeqIO.parse(reads_path, reads_path[-5:]):
        lengths.append(len(record))
    counter = Counter(lengths)
    counter = counter.items()
    counter = sorted(counter)
    with open(profile_path, 'w') as f:
        for pair in counter:
            length, count = pair
            string = str(length) + ', ' + str(count) + '\n' 
            f.write(string)
    print('Generated reads length distribution profile file at ' + os.path.abspath(profile_path))
    return profile_path



def seqreq_perfect_reads(config):
    program_path = config['seqrequester']['program_path']
    genome_path = config['seqrequester']['genome_path']
    coverage = config['seqrequester']['coverage']
    reverse_prob = config['seqrequester']['reverse_prob']

    circular = config['seqrequester'].getboolean('circular')
    length_profile_source = config['seqrequester']['length_profile_source']
    output_path = config['seqrequester']['output_path']
    command = [program_path, 'simulate', '-genome', genome_path, '-coverage', coverage, '-reverse', reverse_prob]
    if circular:
        command.append('-circular')
    command.append('-distribution')
    command.append(get_length_profile(length_profile_source))
    command.append('-genomesize')
    command.append(str(get_genome_length(genome_path)))
    print('Simulating perfect reads...')
    if output_path == '':
        subprocess.run(command)
    else:
        out = open(output_path, 'w')
        subprocess.run(command, stdout=out)
        out.close()
    




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate perfect reads from a genome fasta file.')
    parser.add_argument('-c', '--config', type=str, help='path to the config file.')
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config)
    seqreq_perfect_reads(config)