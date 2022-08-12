# BioDataProcessing
This repo contains scripts to simulate reads with simple substitution, insertion and deletion errors from a genome. 

## Overview
1. simulate_perfect_reads.py - basically just a wrapper over Seqrequester. Requires a config file.

2. add_errors.py - add errors to sequences in a fasta/q file.

3. simulate_reads.py - uses the other scripts to make both reads without errors and with errors given a genome.

------------------------------------------
More details are in the comments and -h.


## Dependencies
### Python libraries
-BioPython

-Numpy

### Tools  
-Seqrequester

