import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import configparser


def rle(seq):
    k, _ = zip(*itertools.groupby(seq))
    return ''.join(k)


def process_fastaq_file(in_path, out_path):
    rles = []
    file_type = in_path[-5:]
    for record in SeqIO.parse(in_path, file_type):
        r = rle(str(record.seq))
        sr = SeqRecord(Seq(r))
        rles.append(sr)

    SeqIO.write(rles, out_path, 'fasta')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run-length encode a FASTA/Q file and output a FASTA file')
    parser.add_argument('-i', '--input', type=str, help='path of the input FASTA/Q file.')
    parser.add_argument('-o', '--output', type=str, help='path of the output FASTA file.')
    args = parser.parse_args()

    process_fastaq_file(args.input, args.output)