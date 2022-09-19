#!/usr/bin/env python

import argparse
import pysam

'''
Inputs: 
1. path to assembly to reference bam file with index.
2. Region, in a tuple (contig, start, end) corresponding to pysam fetch format
3. bool, whether to include the assembly regions that are soft clipped off and consider them as 
mapping to the regions on the ref as if they are all Ms on cigar.
Outputs:
1. Regions of the assembly corresponding to the specified regions on the reference.
In format of a list of (contig, start, end) tuples 
'''

CIGAR_INT2CHAR = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B']
CIGAR_CHAR2INT = {
    'M' : 0,
    'I' : 1,
    'D' : 2,
    'N' : 3,
    'S' : 4,
    'H' : 5,
    'p' : 6,
    '=' : 7,
    'X' : 8,
    'B' : 9,
    'NM' : 10
}
CODE_M = 0
CODE_I = 1
CODE_D = 2
CODE_S = 4
CODE_H = 5

def corresponding_regions(bam_path, region, include_clipped, get_supp = False, get_supp_from_SA=False, supp_include_clipped=False):
    output = []
    
    bam = pysam.AlignmentFile(bam_path) 
    for record in bam.fetch(region[0]):
        if region[2] < record.reference_start or region[1] >= record.reference_end:
            continue
        if record.is_secondary:
            continue

        if record.is_supplementary:
            if not get_supp:
                continue
            #print('is supplementary alignment!')
            segment = infer_corresponding_region(record.query_name, record.reference_start, record.reference_end,
                record.is_reverse, record.cigartuples, region, supp_include_clipped)
        else:
            segment = infer_corresponding_region(record.query_name, record.reference_start, record.reference_end,
                record.is_reverse, record.cigartuples, region, include_clipped)
        output.append(segment)
        if get_supp and get_supp_from_SA and record.has_tag('SA'): # recover filtered supp alignments from SA tag
            supp_alignment_tups = parse_supp_string(record.get_tag('SA'))
            for _, ref_start, ref_end, is_rev, cigar_tuples, _, _ in supp_alignment_tups:
                if region[2] < ref_start or region[1] >= ref_end:
                    continue
                segment = infer_corresponding_region(record.query_name, ref_start, ref_end,
                    is_rev, cigar_tuples, region, supp_include_clipped)
                output.append(segment)
    return output

def find_qpos_change(rpos_change, cigar_tuples):
    matching = 0
    qpos_change = rpos_change

    for c, l in cigar_tuples:
        if c == CODE_M:
            matching += l
            if matching >= rpos_change:
                break
        elif c == CODE_I:
            qpos_change += l
        elif c == CODE_D:
            qpos_change -= min(l, rpos_change-matching)
            matching += l
            if matching >= rpos_change:
                break
        elif c == CODE_S or c == CODE_H:
            continue
        else:
            print('???!')
    return qpos_change
     
def parse_supp_string(string):
    ls = string.split(';')
    ls.pop()
    ls = [x.split(',') for x in ls]
    def parse(x):
        is_rev = False
        if x[2] == '-':
            is_rev = True
        cigar_tuples = []
        start = 0
        current = 0
        ref_start = int(x[1]) - 1
        ref_end = ref_start
        for c in x[3]:
            if c.isalpha():
                op_len = int(x[3][start:current])
                op = x[3][current]
                op = CIGAR_CHAR2INT[op]
                cigar_tuples.append((op, op_len))
                start = current + 1
                if op != CODE_S and op != CODE_I:
                    ref_end += op_len
            current += 1
        return x[0], ref_start, ref_end, is_rev, cigar_tuples, int(x[4]), int(x[5])
    ls = [parse(x) for x in ls]
    return ls


def infer_query_info(cigar_pairs):
    q_len = 0
    seen_M = False
    left_clip = 0
    right_clip = 0
    for op, num in cigar_pairs:
        if op == CODE_M:
            q_len += num
            seen_M = True
        elif op == CODE_I:
            q_len += num
        elif op == CODE_S or op == CODE_H:
            q_len += num
            if seen_M:
                right_clip += num
            else:
                left_clip += num
    
    q_end = q_len - right_clip
    q_start = left_clip

    return q_start, q_end, q_len
def infer_corresponding_region(q_name, ref_start, ref_end, is_rev, cigar_tuples, region, include_clipped):
    q_start, q_end, q_len = infer_query_info(cigar_tuples)
    start = None
    end = None
    if ref_start < region[1]:
        start = q_start + find_qpos_change(region[1]-ref_start, cigar_tuples)
    else:
        start = q_start
        if include_clipped:
            num_left_soft_clip = 0
            if cigar_tuples[0][0] == CODE_S :
                num_left_soft_clip = cigar_tuples[0][1]
            if len(cigar_tuples) > 1 and cigar_tuples[1][0] == CODE_S: # maybe the first one can be hard-clip?..
                num_left_soft_clip = cigar_tuples[1][1]
            start -= min(num_left_soft_clip, ref_start - region[1])
    if ref_end > region[2]:
        end = q_end - find_qpos_change(ref_end - region[2], reversed(cigar_tuples))
    else:
        end = q_end
        if include_clipped:
            num_right_soft_clip = 0
            if cigar_tuples[-1][0] == CODE_S:
                num_right_soft_clip = cigar_tuples[-1][1]
            if len(cigar_tuples) > 1 and cigar_tuples[-2][0] == CODE_S: # maybe the first one can be hard-clip?..
                num_right_soft_clip = cigar_tuples[-2][1]
            end += min(num_right_soft_clip, region[2] - ref_end) 
    if is_rev:
        rstart = q_len - end 
        rend = q_len - start
        return q_name, rstart, rend
    return q_name, start, end





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Deals with regions of alignments.')
    parser.add_argument('-i', '--bam', type=str, help='path of the bam.')
    args = parser.parse_args()

    x = corresponding_regions(args.bam, ('1', 12000, 28000), True)
    print(x)



'''
start = record.query_alignment_start 
end = record.query_alignment_end
print(start, end)
left_offset = record.reference_start - region[1]
right_offset = record.reference_end - region[2]
print(left_offset)
if left_offset < 0: # reference_start excludes clipped
    start -= left_offset
if right_offset > 0:
    end -= right_offset


if include_clipped:
    num_left_soft_clip = 0
    if record.cigartuples[0][0] == 4:
        num_left_soft_clip = record.cigartuples[0][1]
    if len(record.cigartuples) > 1 and record.cigartuples[1][0] == 4: # maybe the first one can be hard-clip?..
        num_left_soft_clip = record.cigartuples[1][1]
    start -= (max(0, min(num_left_soft_clip, left_offset)))            

    num_right_soft_clip = 0
    if record.cigartuples[-1][0] == 4:
        num_right_soft_clip = record.cigartuples[-1][1]
    if len(record.cigartuples) > 1 and record.cigartuples[-2][0] == 4: # maybe the first one can be hard-clip?..
        num_right_soft_clip = record.cigartuples[-2][1]
    end += (max(0, min(num_right_soft_clip, -right_offset)))            
'''