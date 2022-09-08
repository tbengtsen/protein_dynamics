import mdtraj as md
import argparse, math, sys, os, time, glob, json, random 
import numpy as np
import Bio.PDB
# from pdbfixer.pdbfixer import PDBFixer
from dump_json import *
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
                        help='list or json that contains the pdbs names of simulations to dump.  Default ../data/CATH4.2_ALL_STATUS.json ')
    parser.add_argument('-sim_dir',
                        dest='sim_dir',
                        type=str,  
                        required = False,
                        default = '../simulations_new/',
                        help='path to directory with all simulations')
    
    parser.add_argument('-o',
                        dest='out', 
                        type=str,
                        required = False,
                        default='../dynamics_augmented_strucutures/simulated_trajs.jsonl',
                        help='output jsonl file. Default: ../dynamics_augmented_strucutures/simulated_trajs.jsonl')
    
    parser.add_argument('-freq',
                        dest='freq',
                        type=int,
                        default = 1000,
                        required = False,
                        help='how often to output frames from simulation in ps. Default = all' )
    
    parser.add_argument('-max_len',
                        dest='max_len',
                        type=int,
                        default = 500,
                        required = False,
                        help='Filter function. Only output simulations of proteins <=max_len. Default <=500' )
    
    parser.add_argument('-exp_method',
                        dest='exp',
                        type=str,
                        default = 'all',
                        choices = ['Xray','NMR','all'],
                        required = False,
                        help='Filter function. Only output simulations of proteins experimental solved by either Xray crystallization, NMR, or all Default: all' )
    
    parser.add_argument('-max_chains',
                        dest='max_chains',
                        type=int,
                        default = -1,
                        required = False,
                        help='Filter function. Set the maximum number of chains allowed in the protein. Default: all chains (-1)' )
    

    
    
    
    args = parser.parse_args()

    return args


def retrieve_pdb_list(path_list_pdbs,
                     filter_max_len = None,
                     filter_exp_method = None,
                     filter_max_chains = None):
    '''get list of simulated pdbs in the list or json format '''
    
    if path_list_pdbs.split('.')[-1] == 'json':
        with open(path_list_pdbs, 'r') as f:
            all_cath = json.load(f)
        simulated_pdbs = [pdb for pdb in all_cath if all_cath[pdb]['simulation status'] == 'Finished']
        

        
        # filter pdbs if specified
        if  filter_max_len != None:
            simulated_pdbs = [pdb for pdb in simulated_pdbs if all_cath[pdb]['size'] <= filter_max_len]
        
        if filter_exp_method != None and filter_exp_method != 'all' :
            if filter_exp_method == 'Xray': method = 'X-RAY DIFFRACTION'  
            if filter_exp_method == 'NMR' : method = 'SOLUTION NMR' 
            simulated_pdbs = [pdb for pdb in simulated_pdbs if method in all_cath[pdb]['method']]

        if filter_max_chains != -1 and filter_max_chains is not None:
            simulated_pdbs = [pdb for pdb in simulated_pdbs if all_cath[pdb]['chains'] <= filter_max_chains]

            
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

    
    # get list of all simulated pdbs:
    simulated_pdbs = retrieve_pdb_list(args.list_pdbs, 
                                       args.max_len,
                                       args.exp,
                                       args.max_chains)
    
    # collect all dumped 
    all_trajs = []
    
    for pdb in simulated_pdbs:
    
        print(f'\n\n Processing {pdb}')
        try:
            path = args.sim_dir
            pdb_dir = os.path.join(path, pdb)
            traj = os.path.join(pdb_dir, f'{pdb}_trajectory.xtc')
            top_pdb = os.path.join(pdb_dir, f'equil_2_{pdb}.pdb')
    #        renumber_by = f'/home/trz846/protein_dynamics/data/pdbs_cleaned_new/{pdb}_clean.pdb'
            dumped_traj = dump_traj(traj_fn = traj, 
                       pdb_fn = top_pdb,
                       freq = args.freq,
                       pdb_name = pdb,
                       )
        
            all_trajs += dumped_traj
        
        # temporary for working on both servers
        except:
#                 path = args.sim_dir + '/simulations_comp_new/' 
                path = args.sim_dir.split('/')[:-1]
                path = '/'.join(path) 
                path = path + '/simulations_nan_cleaned_by_mahers_script/' 
                pdb_dir = os.path.join(path, pdb)
                traj = os.path.join(pdb_dir, f'{pdb}_trajectory.xtc')
                top_pdb = os.path.join(pdb_dir, f'equil_2_{pdb}.pdb')
        #        renumber_by = f'/home/trz846/protein_dynamics/data/pdbs_cleaned_new/{pdb}_clean.pdb'
                dumped_traj = dump_traj(traj_fn = traj, 
                           pdb_fn = top_pdb,
                           freq = args.freq,
                           pdb_name = pdb)
                all_trajs += dumped_traj

    
    
    with open(args.out, 'w') as f:
        for structure in all_trajs:
            f.write(json.dumps(structure)+'\n')
