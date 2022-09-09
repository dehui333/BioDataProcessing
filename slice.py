#!/usr/bin/env python

import argparse
from Bio import SeqIO
from check_assembly import get_diff
from collections import namedtuple
import os.path
from pathlib import Path
import pysam
import subprocess
import sys

'''
Do slicing and modifications of assembly.

* Merging diffs that are close together may be a thing worth having?
'''
'''
A namedtuple representing a difference between the assembly and reference.
type code:
2: del
8: mismatch
1: ins
10: hp del
'''
CODE_MATCH = 7
CODE_DEL = 2
CODE_MIS = 8
CODE_INS = 1
CODE_HP = 10 # HP DEL
CODE_SOFT = 4
CODE_HARD = 5

diff = namedtuple('diff', ['pos', 'type', 'chunk', 'len'])

OP_DEL = 0
OP_INS = 1
OP_REPLACE = 2
OP_SPLIT = 3
operation = namedtuple('operation', ['pos', 'type', 'length', 'data'])


TYPE_COV = 0 
TYPE_INS = 1
TYPE_DEL = 2
TYPE_MIS = 4
TYPE_HP = 5

'''
Gives list of differences between query and ref.
The ins/del will come after the pos.
Mismatch/match/S/H starts from pos.
'''
def contig_ref_diff(query, ref, rstart, cigartuples):
    ls = []
    qpos = 0 # point to next possible match/mismatch/S/H 
    rpos = rstart # point to next possible match/mismatch
    for case, count in cigartuples:
        if case == 7:
            # equal
            qpos += count  
            rpos += count 
            continue
        if case == 2:
            # del
            hp = False
            deleted = ref[rpos:rpos+count]
            left = ref[rpos-3:rpos]
            right = ref[rpos+count:rpos+count+3]
            if right and right[0] == deleted[-1] and right == right[0] * 3:
                hp = True
            elif left and left[-1] == deleted[0] and left == left[0] * 3:
                hp = True
            if hp:
                ls.append(diff(qpos-1, CODE_HP, deleted, count))
            else:
                ls.append(diff(qpos-1, CODE_DEL, deleted, count))
            rpos += count
            continue
        if case == 8:
            # mismatch
            ls.append(diff(qpos, CODE_MIS, ref[rpos:rpos+count], count))
            qpos += count
            rpos += count
            continue
        if case == 1:
            # ins
            ls.append(diff(qpos-1, CODE_INS, query[qpos:qpos+count], count))
            qpos += count
            continue
        if case == 4:
            # S 
            ls.append(diff(qpos, CODE_SOFT, None, count))
            qpos += count
            continue 
        if case == 5:
            # H
            ls.append(diff(qpos, CODE_HARD, None, count))
            qpos += count
            continue 
    return ls
    

'''
Returns a list of differences between assembly and reference.
Unmapped, secondary and supplementary alignments are skipped.
'''
def assm_ref_diff(assm2ref_bam_path, ref_path):
    samfile = pysam.AlignmentFile(assm2ref_bam_path, 'r')
    ref_contigs = SeqIO.index(ref_path, 'fasta')
    output_dict = {} # dict mapping contig name to differences
    #assm_contigs = SeqIO.index(assm_path, 'fasta')
    for record in samfile.fetch():
        if record.is_unmapped:
            continue
        if record.is_secondary:
            continue
        if record.is_supplementary:
            continue
        #output_dict[record.query_name] = 
        output_dict[record.query_name] = (record.is_reverse, contig_ref_diff(record.query_sequence, ref_contigs[record.reference_name].seq,
            record.reference_start, record.cigartuples))
    return output_dict

'''
for assm2ref
Input: A list of diff namedtuples
Output: A list of operations to modify the contig
'''
def make_operations(diffs, is_reverse, contig_len, skip_mismatch=False, fix_HP=False):
    ls = []
    if is_reverse:
        for diff in diffs:
            pos = contig_len - 1 - diff.pos
            # recheck all cases
            if diff.type == CODE_HP:
                if fix_HP:
                    ls.append(operation(pos-1, OP_INS, diff.len, diff.chunk.reverse_complement()))
                else:
                    ls.append(operation(pos, OP_SPLIT, None, None))    
            elif diff.type == CODE_MIS:
                if not skip_mismatch:
                    ls.append(operation(pos - diff.len + 1, OP_DEL, diff.len, None))
                    ls.append(operation(pos + 1, OP_SPLIT, None, None))
            elif diff.type == CODE_DEL:
                ls.append(operation(pos, OP_SPLIT, None, None))
            elif  diff.type == CODE_INS:
                ls.append(operation(pos-diff.len, OP_DEL, diff.len, None))
                ls.append(operation(pos, OP_SPLIT, None, None))
                
            elif diff.type == CODE_SOFT or diff.type == CODE_HARD:
                ls.append(operation(pos-diff.len+1, OP_DEL, diff.len, None))
            else:
                print('STH WRONG')
    else:
        for diff in diffs:
            if diff.type == CODE_HP:
                if fix_HP:
                    ls.append(operation(diff.pos, OP_INS, diff.len, diff.chunk))
                else:
                    ls.append(operation(diff.pos+1, OP_SPLIT, None, None))    
            elif diff.type == CODE_MIS:
                if not skip_mismatch:
                    ls.append(operation(diff.pos, OP_DEL, diff.len, None))
                    ls.append(operation(diff.pos+diff.len, OP_SPLIT, None, None))
            elif diff.type == CODE_DEL:
                ls.append(operation(diff.pos+1, OP_SPLIT, None, None))
            elif  diff.type == CODE_INS:
                ls.append(operation(diff.pos+1, OP_DEL, diff.len, None))
                ls.append(operation(diff.pos + diff.len+1, OP_SPLIT, None, None))
                
            elif diff.type == CODE_SOFT or diff.type == CODE_HARD:
                ls.append(operation(diff.pos, OP_DEL, diff.len, None))
            else:
                print('STH WRONG')
    return ls

# for reads2assm
def make_operations2(diffs):
    # 0 - pos, 1 - type, 2 - len
    ls = []
    for diff in diffs:
        if diff[1] == TYPE_INS:
            ls.append(operation(diff[0]+1, OP_DEL, diff[2], None))
            ls.append(operation(diff[0] + diff[2]+1, OP_SPLIT, None, None))
        elif diff[1] == TYPE_DEL:
            ls.append(operation(diff[0]+1, OP_SPLIT, None, None))
        elif diff[1] == TYPE_MIS:
            ls.append(operation(diff[0], OP_DEL, diff[2], None))
            ls.append(operation(diff[0]+diff[2], OP_SPLIT, None, None))
        elif diff[1] == TYPE_COV:
            ls.append(operation(diff[0], OP_DEL, diff[2], None))
            ls.append(operation(diff[0]+diff[2], OP_SPLIT, None, None))
        elif diff[1] == TYPE_HP:
            continue # do nothing
        else:
            print('STH WRONG2') 
    return ls
'''
'modifications' is a list of operations.
There are 4 types of operations - deletion, insert, replace and split.
Each operation is represented by a namedtuple operation(pos, type, length, data)
del will delete from pos onwards, ins insert after pos, replace will replace starting from pos,
split will make pos the start of the next chunk.
'''     
def modify_contig(modifications, contig_name, contig):
    chunk_idx = 0
    new_record = True
    pos = 0 # pos idx of next to be output
    # requires a non empty sequence segment as input
    def print_chunk(content, end_of_record=False):
        nonlocal chunk_idx
        nonlocal new_record
        nonlocal contig_name
        if new_record:
            print('>'+contig_name + '_' + str(chunk_idx))
            chunk_idx += 1
            new_record = False
        if end_of_record:
            print(content)
            new_record = True
        else:
            print(content, end='')

    for op in modifications:
        if op.type == OP_INS:
            if pos <= op.pos: 
            # print till include op.pos as ins is after
                print_chunk(contig[pos:op.pos+1])
            print_chunk(op.data)
            pos = op.pos + 1
        elif op.type == OP_REPLACE:
            if pos < op.pos:
                print_chunk(contig[pos:op.pos])
            print_chunk(op.data)
            pos = op.pos + op.length
        elif op.type == OP_DEL:
            if pos < op.pos:
                print_chunk(contig[pos:op.pos])
            pos = op.pos + op.length
        elif op.type == OP_SPLIT: # start new chunk

            if pos < op.pos:
                print_chunk(contig[pos:op.pos], True)
            else:
                if not new_record: # if in the middle of some record
                    print()
                    new_record = True
            pos = op.pos
    if pos < len(contig):
        print_chunk(contig[pos:], True)
    elif not new_record:
        print()



'''
Modify assembly according to its differences with the reference.
contigs: dictionary of fasta records in assembly
diffs: dictionary of contig_name : (is_reverse, list of diffs) 
'''
def modify_assembly(contigs, diffs_assm2ref, r2assm, skip_mismatch=False, fix_HP=False):
    #records =  SeqIO.index(assm_path, 'fasta')
    for contig_name, record in contigs.items():
        ops1 = [] # from assm2ref 
        ops2 = [] # from reads2assm      
        if diffs_assm2ref != None and contig_name in diffs_assm2ref.keys():
            diff_tuple = diffs_assm2ref[contig_name]
            is_reverse = diff_tuple[0]
            differences = diff_tuple[1]
            ops1 = make_operations(differences, is_reverse, len(record), skip_mismatch, fix_HP)
        if r2assm != None:
            diffs = get_diff(r2assm, contig_name, str(record.seq), len(record), 20, 5, 0.4, True)
            op2 = make_operations2(diffs)
        ops = ops1 + ops2
        ops.sort(key=lambda tup: tup[0])
        modify_contig(ops, contig_name, record.seq)     
        
    '''
    for contig_name, diff_tuple in diffs_assm2ref.items():
        #print('modifying ' + contig_name, file=sys.stderr)
        record = contigs[contig_name]
        #print('len: ' + str(len(record)), file=sys.stderr)
        is_reverse = diff_tuple[0]
        differences = diff_tuple[1]
        #if is_reverse:
        #    record = record.reverse_complement()
        ops = make_operations(differences, is_reverse, len(record), skip_mismatch, fix_HP)
        modify_contig(ops, contig_name, record.seq)
    '''

'''
returns sam_path
'''
def align_assm2ref(assm_path, ref_path, num_threads, force):
    assm_dir = os.path.dirname(assm_path)
    if assm_dir != '':
        assm_dir += '/' 
    sam_path =  assm_dir + Path(assm_path).stem + '2' + Path(ref_path).stem + '.sam'
    if force:
        subprocess.run(['rm', sam_path])
    already_has_sam = os.path.exists(sam_path)
    if not already_has_sam:
        with open(sam_path, 'w') as sam_file:
            subprocess.run(['minimap2', '-x', 'asm5', '--secondary=no', '-L', '-a', '--eqx', '-t', \
                str(num_threads), ref_path, assm_path], stdout=sam_file)
    else:
        print(sam_path + ' exists, using.', file=sys.stderr)
    return sam_path




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare assembly with reference/reads.')
    parser.add_argument('-i', '--assm', type=str, help='path of the assembly.')
    parser.add_argument('-r', '--ref', type=str, default=None, help='path of the reference fasta.')
    parser.add_argument('--rd', type=str, default=None, help='Reads to assembly bam.')
    #parser.add_argument('-o', '--output', type=str, help='path of the output assembly fasta.')
    parser.add_argument('-t', '--num_threads', type=str, help='number of threads to use.')
    parser.add_argument('-f', action='store_true', help='force creation of sam.')
    args = parser.parse_args()

    contigs =  SeqIO.index(args.assm, 'fasta')
    diffs_assm2ref = None
    if args.ref != None:
        assm2ref_sam_path = align_assm2ref(args.assm, args.ref, args.num_threads, args.f)
        diffs_assm2ref = assm_ref_diff(assm2ref_sam_path, args.ref)
    modify_assembly(contigs, diffs_assm2ref, args.rd, skip_mismatch=True, fix_HP=True)
