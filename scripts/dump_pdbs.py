import os,glob, json
import numpy as np
import mdtraj as md
import argparse
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
                        help='name of traj file')
    parser.add_argument('-pdb',
                        dest='pdb',
                        type=str,
                        help='path to pdb topology file, (only protein no solvent)'
                       )
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
                        help='name of the outputted pdb frame files (default: <pdb>_XXps.pdb)')
    
    parser.add_argument('-renumber_file',
                        dest='renumber_file', 
                        type=str,
                        help='pdb file to renumber outpdbs after. This is because output files from the openMM simulation does not preserve residue numbering')
    
    args = parser.parse_args()

    return args

    
    
def dump_traj(traj_fn, 
                pdb_fn,
                freq = None,
                o_name=None, 
                o_dir = './traj_pdbs',
                renumber_file = None ):
    '''
    Takes trajectory file and cleaned pdb file and splits
    trajectory into n frames that are each save as a pdb file
    '''
    assert os.path.exists(o_dir), 'Output directory does not exit'
    
    # if top and traj does not match in number of atoms, catch and try to 
    # only load protein in topology
    try:
        traj = md.load(traj_fn, top=pdb_fn)
    except ValueError:
        print('OBS! Number of atoms in trajectory file and pdb file does not',
                'match. Trying to read only protein atoms\n\n')
        all = md.load_pdb(pdb_fn)
        protein = all.topology.select('protein')
        prot_idx = (len(protein))
        top =  md.load_pdb(pdb_fn, atom_indices=np.arange(prot_idx))
        traj = md.load(traj_fn, top=top)

    if o_name is None:
        o_name = pdb_fn.split('/')[-1].split('.pdb')[0]
    

    
    # if no timestep freq, output all frames in traj
    if freq == None: 
        freq = 1 # outputs all frames saved in traj

    for frame in traj:
        ts = int(frame.time[0])
        first_frame = int(traj.time[0])
        last_frame = int(traj.time[-1])
        if ts % freq == 0 or ts==first_frame or ts==last_frame:
            o_pdb = o_dir + f'/{o_name}_{ts}ps.pdb'
            if renumber_file is None :
                frame.save_pdb(o_pdb) 
            else:
                structure = renumber(frame, renumber_file)
                with open(o_pdb, "w") as out:
                    PDBIO.set_structure(structure)
                    PDBIO.save(out)
    
    # also dump pdb database original structure
    o_pdb = o_dir + f'/{o_name}_-1ps.pdb'
    if renumber_file is None :
        top.save_pdb(o_pdb) 
    else:
        structure = renumber(top, renumber_file)
        with open(o_pdb, "w") as out:
            PDBIO.set_structure(structure)
            PDBIO.save(out)



def renumber(frame, renumber_file):
    '''
    Correct for pdbfixer not preserving insertion codes
    '''
    # create temp file to save frame in 
    with tempfile.NamedTemporaryFile(mode="wt", delete=True) as temp:
        frame.save_pdb(temp.name) 
        temp.flush()
        # Fix IDs manually since pdbfixer does not preserve insertion codes
        structure_before = PDB_PARSER.get_structure(renumber_file,renumber_file)
        structure_after = PDB_PARSER.get_structure(temp.name, temp.name)
    residues_before = []
    for chain in structure_before[0]:
        residues_before.append(chain.get_list())
    residues_after = []
    for chain in structure_after[0]:
        residues_after.append(chain.get_list())
    chain_counter = ""
    for i, chain in enumerate(structure_before[0]):
        try:
            if (
                structure_after[0].get_list()[i].id
                != structure_before[0].get_list()[i].id
            ):
                try:
                    # HACK BECAUSE OF https://github.com/biopython/biopython/issues/1551
                    # Essentially, a new change in biopython prevents you from changing the
                    # id to an already existing id which broke this initial script.
                    # Therefore, we now change the ids to "change_counter" which will never look
                    # like a canonical chainid.
                    structure_after[0][
                        structure_before[0].get_list()[i].id
                    ].id = chain_counter
                    chain_counter += "KK"
                except KeyError:
                    pass
                structure_after[0].get_list()[i].id = (
                    structure_before[0].get_list()[i].id
                )
            if len(residues_before[i]) != len(residues_after[i]):
                print('Error in lenght before and after')
                print('before:',len(residues_before[i]), 'after',len(residues_after[i]) )
                raise PDBFixerResIdentifiabilityIssue()

        # When exceeding chainid Z, pdbfixer has discarded it, whereas biopython has not.
        # For simplicity, we just discard it as well and pretend it does not exist.
        # This is a very rare instance and will likely never be a problem unless you
        # are extremely unlucky to work with huge proteins where you care about the
        # truncation.
        except IndexError:
            continue

        counter = 99999  # A large residue number that will never exist in a pdb.
        for res1, res2 in zip(residues_before[i], residues_after[i]):
            assert (
                res1.get_resname().strip() == res2.get_resname().strip()
                or pdbfixer.pdbfixer.substitutions[res1.get_resname()].strip()
                == res2.get_resname().strip()
            )
            if res2.id != res1.id:
                try:
                    # Similar issue as previous hack https://github.com/biopython/biopython/issues/1551
                    structure_after[0][chain.get_id()][res1.id].id = (
                        " ",
                        counter,
                        " ",
                    )
                except KeyError:
                    pass
                res2.id = res1.id
                counter += 1

    return structure_after
        

    
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
        
    dump_traj(traj_fn = args.traj, 
               pdb_fn = args.pdb,
               freq = args.freq,
               o_name = args.o_name,
               o_dir = args.o_dir,
               renumber_file = args.renumber_file)
    
    
    
