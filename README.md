# BioDataProcessing
This repo contains scripts to deal with biological sequence data.

## Overview
1. simulate_perfect_reads.py - basically just a wrapper over Seqrequester. Requires a config file.

2. add_errors.py - add errors to sequences in a fasta/q file.

3. simulate_reads.py - uses the other scripts to make both reads without errors and with errors given a genome.

4. estimate_errors.py - estimate error rates in a set of reads with respect to a reference assembly.
------------------------------------------
More details are in the comments and -h.

## Dependencies
See the top of each script respectively.

