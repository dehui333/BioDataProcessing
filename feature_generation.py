#!/usr/bin/env python

import align_reads_gen
from itertools import chain
import numpy as np
from typing import Sequence

class FeatureGenerator:
    def __init__(self, reads_path:str, truth_path:str, number_of_threads:int, num_steps:int, debug:bool):
        '''
        Initializes the feature generator.

        Parameters:
        reads_path (str) : path to the input reads file.
        truth_path (str) : path to the ground-truth reads file. None if not producing truth matrices.
        number_of_threads (int) : number of threads to use at the C++ level.
        num_steps (int) : number of positions included in one matrix.
        debug (bool) : whether to prepare for debug printing.

        '''
        # passing int to python level as bool
        debug_int = 0
        if debug:
            debug_int = 1

        # using empty string to indicate no truth
        if truth_path == None:
            truth_path = ''
        
        self.num_steps = num_steps
        self.debug = debug
        align_reads_gen.initialize([[reads_path]], truth_path, number_of_threads, num_steps, debug_int)
    
    def has_next(self):
        return align_reads_gen.has_next()
                
    def get_matrices(self, target_number:int): 
        """
        Try to return 'target_number' of feature matrices of each type:

        The actual number returned will likely be higher because the number of samples produced from each request will vary.
        It may also be lower if there are no more target reads.
        """ 
        count = 0
        matrices_tuple_list = []
        while align_reads_gen.has_next() and count < target_number:
            matrices_tuple = align_reads_gen.generate_features()
            matrices_tuple_list.append(matrices_tuple)
            count += len(matrices_tuple[0])

        output_lists = []
        for i in range(len(matrices_tuple_list[0])):
            output_lists.append(list(chain(*[tuple[i] for tuple in matrices_tuple_list])))
        return output_lists
        
    

def concatenate_and_transpose_train(matrices: Sequence[list[np.ndarray]]) -> tuple[np.ndarray, np.ndarray]:
    '''
    Rearranges the matrices for one target read to form:
    1. A matrice containing all input features for one target read, dimension (number of base positions, 10)
    2. A (number of base positions, ) 1d-matrix for the ground-truth
    '''
    counts_with_stats = [np.transpose(np.concatenate((counts, stats), axis=0)) for counts, stats in zip(matrices[1], matrices[2])]
    all_x_concat = np.concatenate(counts_with_stats, axis=0)
    all_y_concat = np.squeeze(np.concatenate(matrices[0], axis=1))
    return all_x_concat, all_y_concat
    
def main():
    gen = FeatureGenerator('/home/lindehui/projects/reads_correction/dataset/packaged/yeast/reads/FSY1742_tagged_clipped_50k.fastq',
                    '/home/lindehui/projects/reads_correction/dataset/packaged/yeast/reads/FSY1742_truth.fasta',
                    30, 100, True)
    while (gen.has_next()):
        x = gen.get_matrices(100)

if __name__ == '__main__':
    main()