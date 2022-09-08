'''
Cleans a downloaded PDB file to be able to parse to MD protocol. 
This includes a multitude of hacks that each decrease the number of time
the parsing crashes when running on a large number of pdb files.
Not all steps are necessary in most cases, and often only the 
Non heterogenselector, the pdbfixer and the last reduce step is necessary. 
Build on Maher Kassem's Cavity_Model_Demo clean_pdb.py but added 
three more steps, including a second and last round of reduce 
to avoid hydrogens with weird positions in aromatics rings, a centering of 
the structure at the origin and a step to bypass problems with rebuilding
gaps in the structure using pdbfixer due to 1st reduce step removing the 
SEQRES of the pdb header causing the pdb fixer not to locate gaps in 
the structure. 
. 
https://github.com/mahermkassem/Cavity_Model_Demo/
'''
import argparse
import os
import subprocess
import sys
import tempfile
import numpy as np
from io import BytesIO, StringIO

import Bio.PDB
import Bio.PDB.Polypeptide
from Bio.PDB.Polypeptide import PPBuilder
import Bio.SeqIO
import pdbfixer
import simtk
import simtk.openmm
import simtk.openmm.app
import simtk.unit as unit


PDBIO = Bio.PDB.PDBIO()
PDB_PARSER = Bio.PDB.PDBParser(PERMISSIVE=0)


class NonHetSelector(Bio.PDB.Select):
    """ 
    Remove HET atoms and choose first conformation of disordered atoms
    """

    def accept_residue(self, residue):
        norm_res_bool = residue.get_resname() in [
            pdbfixer.pdbfixer.substitutions[key]
            for key in pdbfixer.pdbfixer.substitutions
        ]
        abnorm_res_bool = residue.get_resname() in [
            key for key in pdbfixer.pdbfixer.substitutions
        ]
        return (
            norm_res_bool or abnorm_res_bool
        )  # Accept abnorm since they are converted later

    def accept_atom(self, atom):
        return (
            not atom.is_disordered()
            or atom.get_altloc() == "A"
            or atom.get_altloc() == "1"
        ) and atom.id[0] in ["C", "H", "N", "O", "S", "P"]


class PDBFixerResIdentifiabilityIssue(Exception):
    pass


def _locates_gaps(pdb_input_filename):       
    '''Check pdbfixer locates any missing/unsolved residues in structure in the 
    middle of the sequence. The check is based on structure and not on jumps
    in sequence (the latter can happen when lab cuts out part of seq)
    ''' 
    first_model = PDB_PARSER.get_structure(pdb_input_filename, pdb_input_filename)[0]
    
    # use ppbuilder to locate structural gaps
    ppb = PPBuilder()
    for i, chain in enumerate(first_model):
        pps = ppb.build_peptides(chain)
        if len(pps) > 1:
            # chain breaks in structure (not in sequence) if missing residues in pps => pps>1
            gap =  True
            break

        elif len(pps)==1:
            gap =  False
    
    # check if pdbfixer locates missing gaps
    # pdb fixer only locates a gap if there is a SEQRES 
    # in the pdb header
    missing_residues = _get_missing_res(pdb_input_filename)
  

    if missing_residues and gap: 
        return True
    if not missing_residues and not gap:
        return True
    else:
        return False
        

    
    
def _get_missing_res(pdb_input_filename):
    ''' 
    Retrieves missing residues from the SEQRES record from the pdb. 
    This is needed for pdb fixer to fill missing residues. However the 
    reduce processing below looses the SEQRES record, meaning 
    '''
    pre_fixer = pdbfixer.PDBFixer(pdb_input_filename)
    pre_fixer.findMissingResidues()
    
    # only fix missing residues in gaps not in terminals
    # update missing residues to not include terminals
    missing_residues = {}
    
    

    chains = list(pre_fixer.topology.chains())
    keys = pre_fixer.missingResidues.keys()
    for key in keys:
        chain = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain.residues())):
            pass
        else:
            missing_residues[key]= pre_fixer.missingResidues[key]

    return  missing_residues



def _reduce(
    reduce_executable,
    pdb_input_filename,
    pdbid,
    temp1,
):
    '''
    Add hydrogens using reduce program -- needed to not crash in next step. 
    '''
    # check if files exists
    #assert os.path.isfile(pdb_input_filename),  f'Cannot find {pdb_input_filename}'

    
    reduce_het_dict = os.path.join(os.path.dirname(reduce_executable), "reduce_wwPDB_het_dict.txt")
    assert os.path.isfile(reduce_het_dict), (
        f'Cannot find reduce file {reduce_wwPDB_het_dict.txt}',
        'please download and put in same directory as reduce executable')
    
    # reduce command 
    command = [
        reduce_executable,
        "-BUILD",
        "-DB",
        os.path.join(os.path.dirname(reduce_executable), "reduce_wwPDB_het_dict.txt"),
        "-Quiet",
        pdb_input_filename,
    ]
    error_code = subprocess.Popen(command, stdout=temp1).wait()
    temp1.flush()
    
    # get first model (if NMR) of pdb
    first_model = PDB_PARSER.get_structure(pdbid, temp1.name)[0]

    return first_model

def  _non_het_filter(structure, temp2):
    '''
    Applies non hetero atom filter. Return structure. 
    '''
    PDBIO.set_structure(structure)
    PDBIO.save(temp2, select=NonHetSelector())
    temp2.flush()
    # get first model in structure
    structure = PDB_PARSER.get_structure(temp2.name, temp2.name)[0]
    
    return structure 

def _pdbfixer(first_model, temp3, pdb_input_filename):
    '''
    Fixes possible problems in input pdb structure e.g missing atoms 
    in the structure.
    '''
    # get the missing res from the SEQRES record of the orig
    # pdb, this is lost in the processing from reduce.
    missing_res = _get_missing_res(pdb_input_filename)
    
    # used to fix numbering in step4 if no gaps in structure
    for chain in first_model:
        for res in chain:
            for atom in res:
                atom.set_altloc(" ")
    PDBIO.set_structure(first_model)
    PDBIO.save(temp3)
    temp3.flush()

    # Use PDBFixer to fix common PDB errors
    fixer = pdbfixer.PDBFixer(temp3.name)
    fixer.missingResidues = missing_res # taken from pdb input structure
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # track if any gaps are rebuilt - if so, skip step4
    if fixer.missingResidues.keys(): rebuilt_res = True 
    else: rebuilt_res = False
        
    
    
    return temp3, fixer, rebuilt_res


def _center_at_origin(fixer, tempX):
    '''Function to center the protein at origin'''
 
    # get coordinates of all atom positions
    positions = fixer.positions
    
    # convert to numpy to make subtraction of 3D coordinate matrix
    positions = np.array([pos.value_in_unit(unit.nanometer) for pos in positions])
    
    # get center of geometry - mean x,y,z 
    mean_pdb = np.mean(positions, axis=0, keepdims=True)
    
    # translate to origin
    positions_centered = positions  - mean_pdb
    # convert back to openMM quantity object in same units 
    positions_centered = simtk.unit.quantity.Quantity(positions_centered, unit.nanometer)
    
    
    simtk.openmm.app.PDBFile.writeFile(
    fixer.topology, positions_centered, tempX, keepIds=False
    )
    tempX.flush()
    
    return fixer.topology, positions_centered, tempX


def _fix_numbering(topology, positions, temp4, temp5):
    '''
    Correct for pdbfixer not preserving insertion codes
    '''
    simtk.openmm.app.PDBFile.writeFile(
        topology,  positions, temp5, keepIds=False
    )
    temp5.flush()
    # Fix IDs manually since pdbfixer does not preserve insertion codes
    structure_before = PDB_PARSER.get_structure(temp4.name, temp4.name)
    structure_after = PDB_PARSER.get_structure(temp5.name, temp5.name)
    residues_before = []
    for chain in structure_before[0]:
        print('before', chain.get_list())
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



    
def clean_pdb(pdb_input_filename: str, out_dir: str, reduce_executable: str):
    """
    Function to clean pdbs using reduce and pdbfixer. The output file will 
    have "_cleaned.pdb" suffix.
    
    Parameters
    ----------
    pdb_input_filename: str
        PDB filename
    out_dir: str
        Output directory.
    reduce_executable: str
        Path to the reduce executable
    """
    
    #assert _locates_gaps(pdb_input_filename), ('PDBfixer has not located',
    #                                   'a possible gap in the structure.',
    #                                   'This is likely due to a problem',
    #                                   'reading the SEQRES in the pdbheader')

    
    pdbid = pdb_input_filename.split("/")[-1].split(".pdb")[0]

    # Step 1: Add hydrogens using reduce program -- needed to not crash in step2 
    with tempfile.NamedTemporaryFile(mode="wt", delete=True) as temp1:
        structure = _reduce(reduce_executable, pdb_input_filename, pdbid, temp1)

    # Step 2: NonHetSelector filter
    with tempfile.NamedTemporaryFile(mode="wt", delete=True) as temp2:
        PDBIO.set_structure(structure)
        PDBIO.save(temp2, select=NonHetSelector())
        temp2.flush()
        structure = PDB_PARSER.get_structure(temp2.name, temp2.name)[0]

            
    # Step 3: Use pdbfixer to correct errors 
    # in pdb structure and replace altloc chars to " " 
    with tempfile.NamedTemporaryFile(mode="wt", delete=True) as temp3:
        temp_3, fixer, rebuilt_residues = _pdbfixer(structure, 
                                                    temp3,
                                                    pdb_input_filename
                                                   )
               
        ## TONE
        # Step 4: center protein at origin 
        with tempfile.NamedTemporaryFile(mode="wt", delete=True) as temp4:
            topology, positions, tempX = _center_at_origin(fixer,temp4)


            if rebuilt_residues is False: 
                # Step 5: Correct for pdbfixer not preserving insertion codes
                with tempfile.NamedTemporaryFile(mode="wt", delete=True) as temp5:
                    structure = _fix_numbering(topology, positions, temp3, temp5)
                with tempfile.NamedTemporaryFile(mode="wt", delete=True) as temp6:
                    PDBIO.set_structure(structure)
                    PDBIO.save(temp6)
                    temp6.flush()

                    # step 6: Add hydrogens again for the missing residues and 
                    # atoms pdbfixer rebuilt 
                    with tempfile.NamedTemporaryFile(mode="wt", delete=True) as temp7:
                        structure = _reduce(reduce_executable, temp6.name, pdbid, temp7)

            elif rebuilt_residues is True:
                print('true')
                with tempfile.NamedTemporaryFile(mode="wt", delete=True) as temp5:
                    simtk.openmm.app.PDBFile.writeFile(
                        topology, positions, temp5, keepIds=True)
                    temp5.flush()
                    # step 6: Add hydrogens again for the missing residues and 
                    # atoms pdbfixer rebuilt 
                    with tempfile.NamedTemporaryFile(mode="wt", delete=True) as temp6:
                        structure = _reduce(reduce_executable, temp5.name, pdbid, temp6)


    with open(
        os.path.join(out_dir, pdbid + "_clean.pdb"), "w"
    ) as outpdb:
        PDBIO.set_structure(structure)
        PDBIO.save(outpdb)


if __name__ == "__main__":
    # Argument Parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_file_in", type=str)
    parser.add_argument("--out_dir", type=str)
    parser.add_argument("--reduce_exe", type=str)

    # Parse arguments
    args_dict = vars(parser.parse_args())
    pdb_input_filename = args_dict["pdb_file_in"]
    out_dir = args_dict["out_dir"]
    reduce_executable = args_dict["reduce_exe"]

    # Clean 
    clean_pdb(pdb_input_filename, out_dir, reduce_executable)
