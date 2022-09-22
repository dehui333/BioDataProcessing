#!/usr/bin/env python
import os
import subprocess


'''
Helper functions to help run programs with the python subprocess module.
program - program filename. Expected to be in PATH.
arguments - a dictionary where the key value pairs define the arguments and options.
'''

def run(program, paired_arguments, single_arguments, redirect_stdout=None):
    ls = [program]
    if paired_arguments != None:
        for key, value in paired_arguments.items():
            ls.append(key)
            ls.append(value)
    
    if single_arguments != None:
        for argument in single_arguments:
            ls.append(argument)
    if redirect_stdout:
        with open(redirect_stdout, 'w') as handle:
            subprocess.run(ls, stdout=handle)
    else:
        subprocess.run(ls)

'''
Template:
def run_xx():
    paired_arguments = {
        '-o' : output,
    }

    single_arguments = [
        in1,
        in2
    ]

    run(xx, paired_arguments, single_arguments)
'''

def run_hifiasm_hetero_reads_only(output_prefix, num_threads, list_of_reads_paths, reuse):
    hap1_fasta = output_prefix + '.hap1.fasta'
    hap2_fasta = output_prefix + '.hap2.fasta'
    if reuse and os.path.isfile(hap1_fasta) and os.path.isfile(hap2_fasta):
        return hap1_fasta, hap2_fasta
    named_arguments = {
        '-o' : output_prefix,
        '-t' : str(num_threads)
    }

    run('hifiasm', named_arguments, list_of_reads_paths)
    hap1_gfa = output_prefix + '.bp.hap1.p_ctg.gfa'
    hap2_gfa = output_prefix + '.bp.hap2.p_ctg.gfa'
    # awk '/^S/{print ">"$2;print $3}' test.p_ctg.gfa > test.p_ctg.fa
    run('awk', None, ['/^S/{print ">"$2;print $3}', hap1_gfa], hap1_fasta)
    run('awk', None, ['/^S/{print ">"$2;print $3}', hap2_gfa], hap2_fasta)
    return hap1_fasta, hap2_fasta

def run_pomoxis_assess_assm(assm_path, ref_path, num_threads, output_prefix='assm'):
    named_arguments = {
        '-i' : assm_path,
        '-r' : ref_path,
        '-d' : 'map-ont',
        '-c' : str(0),
        '-a' : '',
        '-t' : str(num_threads),
        '-p' : output_prefix
    }
    run('assess_assembly', named_arguments, None)

def run_quast(output_dir, assm_paths, ref_path, num_threads):
    named_arguments = {
        '-o' : output_dir,
        '-r' : ref_path,
        '-t' : str(num_threads)
    }

    run('quast.py', named_arguments, assm_paths)

if __name__ == '__main__':
    run_pomoxis_assess_assm('test.asm.hap1.fasta', 'data/references/S288C_reference.fasta', 20)
