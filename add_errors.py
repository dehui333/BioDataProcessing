import argparse
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from functools import partial
from multiprocessing import Pool
import numpy as np
import random

'''
This is a script to mutate (add errors) to sequences in fasta/q files.
I would have written it in C++ except that that's too troublesome.
I hope this is not too slow.

Dehui 10/08/2022
'''

AlPHABET = ['A', 'C', 'G', 'T']
BASE_TO_INT = [
None, None, None, None, None, None, None, None, None, None,
None, None, None, None, None, None, None, None, None, None,
None, None, None, None, None, None, None, None, None, None,
None, None, None, None, None, None, None, None, None, None,
None, None, None, None, None, None, None, None, None, None,
None, None, None, None, None, None, None, None, None, None,
None, None, None, None, None, 0, None, 1, None, None,
None, 2, None, None, None, None, None, None, None, None,
None, None, None, None, 3 
]


'''
Returns 4 lists. The first 2 are of the same lengths.
1. The error positions, 
- 0 indexed
2. The error types,
- 0 for substitution, 1 for insertion, 2 for deletion
3. The stream of random bases for insertion
4. The number of shifts(1-3) to transform bases for substitution, according to a
cyclic order. i.e. A -> C -> G -> T -> A -> ...

'''
def generate_errors(seq_len, sub_prob, ins_prob, del_prob):
    # generate a total error count
    total_prob = sub_prob + ins_prob + del_prob
    mean_error_count = seq_len * total_prob
    total_error_count = np.random.normal(mean_error_count, 0.03 * mean_error_count)
    total_error_count = int(total_error_count)
    if total_error_count <= 0:
        return None

    # generate the error posititions
    error_positions = np.random.choice(seq_len, total_error_count, False)
    error_positions.sort()

    # Generate list of errors
    probs = np.array([sub_prob, ins_prob, del_prob])
    # 0 for subs, 1 for ins, 2 for del
    error_types = np.random.choice(np.arange(3), total_error_count, True, probs/np.sum(probs))
    error_counts = Counter(error_types) 
    # Generate the stream of bases for substitution and insertion
    ins_base_stream = np.random.choice(AlPHABET, error_counts[1], True)

    # substitution error is encoded by the numbers 1 to 3 - representing the shift
    # within the bases A C G T, with wrap around   
    subs_shift_stream = np.random.choice([1, 2, 3], error_counts[0], True)
    return error_positions, error_types, ins_base_stream, subs_shift_stream

'''
Applies the specified errors on a sequence, creating a mutated copy.
'''
def apply_errors(seq, error_positions, error_types, ins_base_stream, subs_shift_stream):
    new_string = ''
    ins_base_index = 0
    subs_shift_index = 0
    next_segment_start = 0
    
    for i in range(len(error_positions)):
        error_position = int(error_positions[i])
        base = seq[error_position]
        new_string += seq[next_segment_start:error_position]
        if error_types[i] == 0: # sub
            new_string += AlPHABET[(BASE_TO_INT[ord(base)] + subs_shift_stream[subs_shift_index]) % 4]
            subs_shift_index += 1
        elif error_types[i] == 1: # ins
            new_string += base + ins_base_stream[ins_base_index]
            ins_base_index += 1
             
        next_segment_start = error_position + 1
    new_string += seq[next_segment_start:]
    return new_string

def add_errors(input_path, output_path, sub_prob, ins_prob, del_prob, num_proc):
    file_type = input_path[-5:]
    global f
    def f(seq):
        simulated_errors = generate_errors(len(seq), sub_prob, ins_prob, del_prob)
        if simulated_errors is None:
            return None
        return apply_errors(seq, *simulated_errors)
    if num_proc > 1:
        pool = Pool(num_proc)

    records =  list(SeqIO.parse(input_path, file_type))
    print('Simulating errors...')
    seqs = [record.seq for record in records]

    # plain map seems to be best
    if num_proc > 1:      
        result = pool.map(f, seqs)
        pool.close()
    else:
        result = list(map(f, seqs))    
    for i in range(len(records)):
        if result[i] == None:
            continue
        records[i].seq = Seq(result[i])
    SeqIO.write(records, output_path, file_type)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add simple errors to reads in a fasta/q file.')
    parser.add_argument('-i', '--input', type=str, help='path of the input fasta/q file.')
    parser.add_argument('-o', '--output', type=str, help='path of the output file.')
    parser.add_argument('-t', '--procs', type=int, default = 1, help='number of processes to use.')
    parser.add_argument('--sub', type=float, default=0.03, help='substitution rate.')
    parser.add_argument('--ins', type=float, default=0.03, help='insertion rate.')
    parser.add_argument('--dele', type=float, default=0.03, help='deletion rate.')
    args = parser.parse_args()

    add_errors(args.input, args.output, args.sub, args.ins, args.dele, args.procs)
    

