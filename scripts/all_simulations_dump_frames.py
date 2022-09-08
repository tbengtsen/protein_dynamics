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



def parser():
    '''
    parser for script
    '''
    
    parser = argparse.ArgumentParser(description='Converts trajectories of pdbs to single pdbs for each defined timestep in traj')
    parser.add_argument('-list_pdbs',
                        dest='list_pdbs',
                        type=str,  
                        required = False,
                        default='../data/CATH4.2_ALL_STATUS.json',
                        help='list or json that contains the pdbs names of simulations to dump. Json format: ["test"]["pdb"]. Default ../data/CATH4.2_ALL_STATUS.json  ')
    parser.add_argument('-sim_dir',
                        dest='sim_dir',
                        type=str,  
                        required = False,
                        default = '../simulations_new/',
                        help='path to directory with all simulations. Default = ../simulations_new/')
    parser.add_argument('-o_dir',
                        dest='o_dir', 
                        type=str,
                        required = False,
                        default='./traj_pdbs',
                        help='output directory (default: traj_pdbs)')
    parser.add_argument('-freq',
                        dest='freq',
                        type=int,
                        default = 50,
                        required = False,
                        help='how often to output frames from simulation in ps. Default = all' )

    
    args = parser.parse_args()

    return args


def retrieve_pdb_list(path_list_pdbs):
    '''get the pdbs in the list or json format '''
    
    
    if path_list_pdbs.split('.')[-1] == 'json':
        with open(path_list_pdbs, 'r') as f:
            all_cath = json.load(f)
        simulated_pdbs = [pdb for pdb in all_cath if all_cath[pdb]['simulation status'] == 'Finished']
    
    else:
        with open(path_list_pdbs, 'r') as f:
            lines = f.readlines()
        simulated_pdbs = [pdb.strip() for pdb in lines ]
        
    assert len(simulated_pdbs)> 0, 'Cannot retrieve any pdbs from input -pdbs argument. Change format please'  

    simulated_pdbs = [pdb.lower() for pdb in simulated_pdbs]

    
    return simulated_pdbs
        
        

    

if __name__=='__main__':
    
    args = parser()
    
    # check input arguments 
    assert os.path.isdir(args.sim_dir), 'Directory with simulations specified by -sim_dir does not exists'
    assert os.path.exists(args.list_pdbs), 'the input list or json given to -list_pdbs does not exists '

    
    simulated_pdbs = retrieve_pdb_list(args.list_pdbs)
    
    print(f'\nDumping {len(simulated_pdbs)} simulations every {args.freq} ps')
    
    # create output directory
    if not os.path.isdir(args.o_dir):
        os.mkdir(args.o_dir)
    print(f'\nDumped pdb frames from simulation to be found in {args.o_dir}\n\n')
    
    for pdb in simulated_pdbs:
        already_exists = glob.glob(f'{args.o_dir}/{pdb}*.pdb')
        if already_exists:
            print(f'{pdb} already dumped - skipping')
        if not already_exists:
            print(f'\n\n Processing {pdb}')
            try:
                pdb_dir = os.path.join(args.sim_dir, pdb)
                traj = os.path.join(pdb_dir, f'{pdb}_trajectory.xtc')
                top_pdb = os.path.join(pdb_dir, f'equil_2_{pdb}.pdb')
            #    renumber_by = f'/home/trz846/protein_dynamics/data/pdbs_cleaned_new/{pdb}_clean.pdb'
                dump_traj(traj_fn = traj, 
                           pdb_fn = top_pdb,
                           freq = args.freq,
                           o_name = pdb,
                           o_dir = args.o_dir)
            #               renumber_file = renumber_by)
            
            # following too be deleted later
            except: 
                path = '/home/trz846/pd/simulations_new/simulations_comp_new/'
                pdb_dir = os.path.join(path, pdb)
                traj = os.path.join(pdb_dir, f'{pdb}_trajectory.xtc')
                top_pdb = os.path.join(pdb_dir, f'equil_2_{pdb}.pdb')
                dump_traj(traj_fn = traj, 
                           pdb_fn = top_pdb,
                           freq = args.freq,
                           o_name = pdb,
                           o_dir = o_dir)

