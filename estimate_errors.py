import argparse
import os.path
from pathlib import Path
import pysam
import subprocess

'''
Dependencies:

Python libraries:
-Pysam

Tools:
-Minimap2

This script is for estimating error rates in reads w.r.t. to a reference.
(Technically an assembly can be used as 'reads' too..)

It seems that as the error rates get larger, we tend to overestimate substitution rate
and underestimate insertion and deletion rates? This has been observed in 
pomoxis scripts' assess_assembly. (could be due to misalignments at higher error rate?)
'''


'''
Given a path to a bam/sam file which contains the reads aligned to the reference, output errors rates.
'''
def estimate_error_bam(bam_path):
    total_len = 0
    total_mismatch = 0
    total_ins = 0
    total_del = 0
    bam = pysam.AlignmentFile(bam_path) 
    for alignment in bam.fetch():
        stats = alignment.get_cigar_stats()
        counts = stats[0] 
        total_mismatch += counts[8]
        total_ins += counts[1]
        total_del += counts[2]
        total_len += alignment.query_alignment_length
    bam.close()
    return total_mismatch/total_len, total_ins/total_len, total_del/total_len

'''
Given a path to a read set (fasta/q) and a reference assembly, estimate error rate in read set.
'''
def estimate_error_reads(reads_path, reference_path, keep_sam, num_threads=1): 
    sam_path = os.path.dirname(reads_path) + '/'+ Path(reads_path).stem + '_to_' + Path(reference_path).stem + '.sam'
    already_has_sam = os.path.exists(sam_path)
    if not already_has_sam:
        with open(sam_path, 'w') as sam_file:
            subprocess.run(['minimap2', '-a', '--eqx', '-t', \
                str(num_threads), reference_path, reads_path], stdout=sam_file, stderr=subprocess.DEVNULL)
        
    rates = estimate_error_bam(sam_path)
    if not keep_sam and not already_has_sam:
        subprocess.run(['rm', sam_path])
    return rates



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Estimate the error rate in a set of reads.')
    parser.add_argument('-i', '--input', type=str, help='path of the input fasta/q file.')
    parser.add_argument('-r', '--ref', type=str, help='path of the reference assembly.')
    parser.add_argument('-t', '--num_threads', type=int, default=1, help='number of threads for mapping.')
    parser.add_argument('-k', action='store_true', help='whether to keep sam.')
    args = parser.parse_args()
    subs_rate, ins_rate, del_rate = estimate_error_reads(args.input, args.ref, args.k, args.num_threads)
    print(f'substitution rate: {subs_rate * 100:.3}%')
    print(f'insertion rate: {ins_rate * 100:.3}%')
    print(f'deletion rate: {del_rate * 100:.3}%')