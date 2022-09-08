#!/usr/bin/env python
'''
Created by Tone Bengtsen July 15 2021.
Renumbers residue ids in pdb for a given chain. 
Used for stability calculations where identical numbering of residues is crucial.
Obs, the dump_pdb or dump_json module does the same by modelling after a pdb
with correct residue id. If such pdb does not exist (i.e does not fit mutational data residue numbering) then use this script to create such pdb to renumber after in the dump modules. 

'''

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import sys, argparse

def parser():
    '''
    makes parser for script
    '''
    
    parser = argparse.ArgumentParser(description='Renumber the residues in a given chain in a pdb file. Used for stability calculations where identical numbering of residues is crucial')
    parser.add_argument('-pdb',
                        dest='pdb',
                        type=str,
                        required = True,
                        help='path to pdb  file that is to be renumbered in residue (amino acids) IDs '
                       )
    
    parser.add_argument('-start_resid',
                        dest='start_resid',
                        type=int,
                        required = True,
                        help='chosen residue id (counter) for first residue in pdb  ' )
    
    parser.add_argument('-chain',
                        dest='chain',
                        type=str,
                        required = False,
                        default = 'A',
                        help='chain in pdb that is to be renumbered ' )
    

    parser.add_argument('-o_name',
                        dest='o_name',
                        type=str,
                        required = False,
                        help='name of renumbered output pdb. Default = <input_file>_renumbered.pdb')
    
    args = parser.parse_args()

    return args


if __name__=='__main__':
    
    # load arguments
    args = parser()

    print(f'\nRenumbering pdb file:  {args.pdb} to start residue at {args.start_resid}')
    
    # create new pdb file to not overwrite
    if args.o_name is None:
        pdb_out  = args.pdb.replace(".pdb","") + "_renumbered.pdb"
    else: 
        pdb_out = args.o_name
        
    # load structure:
    parser = PDBParser()
    structure = parser.get_structure(args.pdb.replace(".pdb", ""), args.pdb)
    model = structure[0] # if NMR structure, load best model i.e. often nr 1. 
    chain = model[args.chain]
    
    
    # fix if new start id is higher that original 
    # this will give an error as it causes the same resid in 
    # the structure twice
    if args.start_resid >  int(chain[1].id[1]):
        # fix by tmp artificial high resid number
        tmp_resid = 500 # no proteins alowed larger than 500 here
        for residue in chain:
            residue.id = (' ', tmp_resid , ' ')
            tmp_resid  += 1
    
    
    # renumber residue IDs
    start_resid = args.start_resid
    for residue in chain:
        residue.id = (' ', start_resid , ' ')
        start_resid  += 1

    
    # save:
    w = PDBIO()
    w.set_structure(structure)
    w.save(pdb_out)
        
    
    print(f'Renumbered pdb file saved at {pdb_out}')

