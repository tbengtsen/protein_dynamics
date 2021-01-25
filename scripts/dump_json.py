import mdtraj as md
import argparse, math, sys, os, time, glob, json
import numpy as np

import Bio.PDB
import tempfile
import pdbfixer
from pdbfixer.pdbfixer import PDBFixer
PDBIO = Bio.PDB.PDBIO()
PDB_PARSER = Bio.PDB.PDBParser(PERMISSIVE=0)



class PDBFixerResIdentifiabilityIssue(Exception):
    pass


def parser():
    '''
    makes parser for script
    '''
    
    parser = argparse.ArgumentParser(description='Converts trajectory to single framed pdbs')
    parser.add_argument('-traj',
                        dest='traj',
                        type=str,  
                        required = True,
                        help='name of traj file')
    parser.add_argument('-pdb',
                        dest='pdb',
                        type=str,
                        required = True,
                        help='path to pdb topology file, (only protein no solvent)'
                       )
    parser.add_argument('-pdb_name',
                        dest='pdb_name',
                        required = True,
                        type=str,
                        help='pdb id')
    
    parser.add_argument('-freq',
                        dest='freq',
                        type=int,
                        help='how often to output frames from simulation in ps.Default = all' )
    parser.add_argument('-o_dir',
                        dest='o_dir', 
                        type=str,
                        default='./traj_pdbs',
                        help='output directory (default: traj_pdbs)')

    parser.add_argument('-o_name',
                        dest='o_name',
                        type=str,
                        help='name of saved putput file in jsonl format')
    
    parser.add_argument('-renumber_file',
                        dest='renumber_file', 
                        type=str,
                        help='pdb file to renumber outpdbs after. This is because output files from the openMM simulation does not preserve residue numbering')
    
    args = parser.parse_args()

    return args

# use to change tri-letter abr from structure to single letter
resname_abbr = {
            'MET':'M',
            'ARG':'R',
            'HIS':'H',
            'LYS':'K',
            'ASP':'D',
            'GLU':'E',
            'SER':'S',
            'THR':'T',
            'ASN':'N',
            'GLN':'Q',
            'CYS':'C',
            'GLY':'G',
            'PRO':'P',
            'ALA':'A',
            'VAL':'V',
            'ILE':'I',
            'PHE':'F',
            'TYR':'Y',
            'TRP':'W',
            'LEU':'L'}




def get_coords(frame):
    '''
    retrives the atomic coordinates of the backbone atoms 
    for each chain in each frame 
    '''
    
    # use to get set alphabet naming of chains
    alphabet = list('abcdefghijklmnopqrstuvwxyz')
    
    # get topology of frame
    topology = frame.topology
    
    # get indeces of all backbone atoms (can't set for specific chain )
    CA_idx = topology.select('protein and name CA')
    N_idx = topology.select('protein and name N')
    C_idx = topology.select('protein and name C')
    O_idx = topology.select('protein and name O')
    
    
    coords = {}
    
    for chain_idx, chain in enumerate(topology.chains):
        # get indices of all atoms in the given chain
        all_atoms_idx = [atom.index for atom in topology.chain(chain_idx).atoms ]
    
        # get indices of backbones atoms in specific chain only
        CA_chain_idx = [ atom_idx for atom_idx in all_atoms_idx if atom_idx in CA_idx ]
        N_chain_idx = [ atom_idx for atom_idx in all_atoms_idx if atom_idx in N_idx ]
        C_chain_idx = [ atom_idx for atom_idx in all_atoms_idx if atom_idx in C_idx ]
        O_chain_idx = [ atom_idx for atom_idx in all_atoms_idx if atom_idx in O_idx ]
      
        # get coords for each frames in traj for each backbone atom
        CA_coords = frame.atom_slice(CA_chain_idx).xyz.tolist()[0]
        N_coords = frame.atom_slice(N_chain_idx).xyz.tolist()[0]
        C_coords = frame.atom_slice(C_chain_idx ).xyz.tolist()[0]
        O_coords = frame.atom_slice(O_chain_idx).xyz.tolist()[0]
        chain_coords = {'N':N_coords,'CA':CA_coords, 'C':C_coords, 'O': O_coords}
        coords[alphabet[chain_idx]] = chain_coords
        
    return coords 



def get_seq (topology):
    '''get the primary sequence of structure
    '''
    alphabet = list('abcdefghijklmnopqrstuvwxyz')

    sequences = {}
    for chain_idx, chain in enumerate(topology.chains):
        seq =  [resname_abbr[res.name] for res in topology.chain(chain_idx).residues ]
        seq = ''.join(seq)
        
        sequences[alphabet[chain_idx]] = seq

    return sequences


        
    
def dump_traj(traj_fn, 
                pdb_fn,
                freq = None,
                pdb_name = None,
                renumber_file = None ):
    '''
    Takes trajectory file and cleaned pdb file and splits
    trajectory into n frames and c chains. 
    E.g. dumps seperately for each chain in structure all 
    n_frames structures - to fit with transformer that only works 
    with single chains. 
    '''
    # save all 
    dump_traj = []
    
    # if top and traj does not match in number of atoms, catch and try to 
    # only load protein in topology
    try:
        traj = md.load(traj_fn, top=pdb_fn)
        
    except ValueError:
        print('OBS! Number of atoms in trajectory file and pdb file does not',
                'match. Trying to read only protein atoms\n\n')
        all_atoms = md.load_pdb(pdb_fn)
        protein = all_atoms.topology.select('protein')
        prot_idx = (len(protein))
        top =  md.load_pdb(pdb_fn, atom_indices=np.arange(prot_idx))
        traj = md.load(traj_fn, top=top)
    
    # get name of pdb 
    if pdb_name is None:
        pdb_name = pdb_fn.split('/')[-1].split('.pdb')[0]
        if len(list(pdb_name)) > 4:
            pdb_name = pdb_name.split('_')[0]
        assert len(list(pdb_name)) == 4, 'pdb id name is not 4 characters. Specify name in arguments'
    
    # if no timestep freq, output all frames in traj
    if freq == None: 
        freq = 1 # outputs all frames saved in traj
    
    # get sequences of each chain
    sequences = get_seq (traj.topology)
    
    # save structures for each chain seperately 
    for chain in sequences.keys():
        
                
        for frame in traj:
            # get sim timestep of each frame
            ts = int(frame.time[0])
            first_frame = int(traj.time[0])
            last_frame = int(traj.time[-1])
            
            # create dict for chain/frame structure info in
            dump = {}
            dump["ts"] = ts
            name = pdb_name + '.' + chain.upper() + '.' + str(ts)
            dump["name"] = name 
            dump["seq"] = sequences[chain]
            
            if ts % freq == 0 or ts==first_frame or ts==last_frame:
                coords = get_coords(frame)[chain]
                assert len(coords['N'])==len(sequences[chain]), 'coordinates and sequence len does not match'
                dump["coords"] = coords
            
                dump_traj.append(dump)
    
    return dump_traj    
        

    
if __name__=='__main__':
    
    args = parser()
    
        
    print(f'\nSplitting trajectory {args.traj}.')
    print(f'Dumping frames every {args.freq} ')
    # create output directory
    if not os.path.isdir(args.o_dir):
        os.mkdir(args.o_dir)
    print(f'\nOutput pdb frames to be found in {args.o_dir}\n\n')
    
    if not args.renumber_file is None:
        assert os.path.exists(args.renumber_file)
        print(f'\nRenumbering frames after pdb {args.renumber_file}\n\n')
        
    dumped_traj = dump_traj(traj_fn = args.traj, 
               pdb_fn = args.pdb,
               freq = args.freq,
               pdb_name = args.pdb_name,
               renumber_file = args.renumber_file)
#     with open(outfile, 'w') as f:
#         for entry in dataset:
#             f.write(json.dumps(entry) + '\n')    
    
    if args.o_name is not None:
        o_file = args.o_name
    else:
        o_file = f"dumped_{args.pdb.split('.pdb')[0]}.jsonl"
        
    with open(f"{args.pdb_name}.json", 'w') as f:
        for structure in dumped_traj:
            f.write(json.dumps(structure)+'\n')
    
