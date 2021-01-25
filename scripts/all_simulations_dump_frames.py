import os,glob, json
import numpy as np
import mdtraj as md
import argparse
import Bio.PDB
import tempfile
from pdbfixer.pdbfixer import PDBFixer
from dump_pdbs import *

PDBIO = Bio.PDB.PDBIO()
PDB_PARSER = Bio.PDB.PDBParser(PERMISSIVE=0)



# obtain the desired pdbs with trajs to split
with open('tmp1','r') as f:
    simulated_pdbs = f.readlines()
simulated_pdbs = [pdb.strip() for pdb in simulated_pdbs]

## settings 
path = '/home/trz846/pd/simulations_new/NMRs/'
o_dir = '/home/trz846/protein_dynamics/dynamics_augmented_strucutures/NMRs/'
if not os.path.exists(o_dir):
    os.mkdir(o_dir)

freq = 1000 # unit = ps=> 1 ns
    
for pdb in simulated_pdbs:
    print(f'\n\n Processing {pdb}')
    pdb_dir = os.path.join(path, pdb)
    traj = os.path.join(pdb_dir, f'{pdb}_trajectory.xtc')
    top_pdb = os.path.join(pdb_dir, f'equil_2_{pdb}.pdb')
#    renumber_by = f'/home/trz846/protein_dynamics/data/pdbs_cleaned_new/{pdb}_clean.pdb'
    dump_traj(traj_fn = traj, 
               pdb_fn = top_pdb,
               freq = freq,
               o_name = pdb,
               o_dir = o_dir)
#               renumber_file = renumber_by)


