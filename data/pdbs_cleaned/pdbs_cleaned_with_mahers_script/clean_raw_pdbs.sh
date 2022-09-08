#!/bin/bash

PYTHON=/home/trz846/../chz526/software/miniconda3/envs/py3.6-ddg/bin/python3.6

# Source environment
source activate /home/trz846/../chz526/software/miniconda3/envs/py3.6-ddg

# Reduce executable
reduce_exe="/home/trz846/pd/scripts/reduce/reduce"

# Clean pdbs
for pdb in `cat RMSDs_to_high.dat`;
do
    echo $pdb
    $PYTHON mahers_clean_pdb.py --pdb_file_in ../pdbs_raw/${pdb}.pdb --reduce_exe $reduce_exe --out_dir .
done
