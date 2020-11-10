'''
Functions used to filter out proteins from the CATH4.3
Proteins are removed if they are unsuitable for simulations based on the below: 
    - membrane protein
    - low resolution >2.5 AA
    - unsuitable experimental method
    - unresolved residues inside the sequence (gaps)
    - too big > 500 residues (too expensive to simulate)
    - too small <50 residues (protein does not have innate structure)
    - any peptide < 50 (avoids peptides without structure by themself )
    - unknown/non-standard amino acid
    - unsuitable ligand (e.g lipid)

'''

# to check if chain is broken/ have gap
import Bio.PDB.mmtf as BioMMTF
from Bio.PDB.Polypeptide import PPBuilder

# download/read pdb in mmtf format 
from collections import defaultdict
sys.path.insert(0, '/home/trz846/graph_transformers/')
from utils.data_mmtf import *


def load_membrane_protein_list():
    '''
    Read membrane protein database into list of pdb IDs.
    Using Stephen White's MPstruc database for finding membrane proteins: 
    https://blanco.biomol.uci.edu/mpstruc/?expandTblOnLoad=1
    load xml over all protein structures that are membrane proteins
    '''
    
    import xml.etree.ElementTree as ET
    mpstruc = 'helper_files/membrane_protein_database_mpstruc.xml'
    tree = ET.parse(mpstruc)
    root = tree.getroot()
    # load pdb IDs of all  membrane proteins
    membrane_proteins = [pdb.text for pdb in root.iter('pdbCode')]
    
    
    return membrane_proteins



def membrane_protein(pdb_name):
    '''check if categorised as membrane protein in MPStruc'''
    
    pdb_name = pdb_name.upper() 
    if pdb_name in membrane_proteins:
        return True
    else:
        return False

def unsuitable_experimental_method(mmtf_record):
    '''check that method used to solve protein structure is allowed''' 
    # define allowed methods == good for MD simulations
    allowed_exp = ['ELECTRON MICROSCOPY', 
                    'SOLID-STATE NMR', 
                    'SOLUTION NMR',
                    'X-RAY DIFFRACTION',
                    'ELECTRON CRYSTALLOGRAPHY']

    # get list of methods used for solving structure
    used_exp = mmtf_record.experimental_methods
    
    # see if any of them is among suitable. 
    is_allowed = [exp for exp  in used_exp if exp in allowed_exp]

    if  is_allowed: return False
    else:  return True 


def low_resolution(mmtf_record):
    '''check if resolution lower than 2.5 AA if experimental method is x-ray crystallization'''
    
    try: # can only get resolution for some exp methods 
        resolution = mmtf_record.resolution 
        if resolution > 2.5: return True
        else: return False
    
    except: # lazy as ... I know... 
        return False 



def polymeric_with_single_chain_too_short(mmtf_record):
    '''returns false if polymeric i.e. multiple chains 
    OBS! does not work for NMR - it reads each model as 1 chain '''
    # if homo-polymer then each of the identical chains are put in 1 entity, i.e. 1 entity contains all homomers in pdb
    protein_entities = [entity for entity in  mmtf_record.entity_list if entity['type'] == 'polymer']
    assert len(protein_entities) > 0, f'Error, pdb has 0 protein chains'
    any_chain_too_short = False
    for chain in protein_entities:
        if len(chain['sequence']) < 50:
            any_chain_too_short = True
    return any_chain_too_short

    


def polymeric(biopython_structure):
    count_chains = len(biopython_structure)
    return count_chains


def aa_selenocystein(mmtf_record):

    protein_entities = [entity for entity in  mmtf_record.entity_list if entity['type'] == 'polymer']
    for chain in protein_entities:
        sequence = chain['sequence']
        seleno_cys_present = False
        if 'U' in sequence:
            seleno_cys_present = True
            break
    
    return seleno_cys_present 


def unknown_aa(mmtf_record):
    
    known_aa = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','Y','V','W', 'U'] # U = Selenocystein.
    bad_aa = False
    protein_entities = [entity for entity in  mmtf_record.entity_list if entity['type'] == 'polymer']
    for chain in protein_entities:
        sequence = chain['sequence']
        bad_aa = [aa for aa in sequence if aa not in known_aa]
        if len(bad_aa) > 0: 
            bad_aa = True
            break 
            
    if bad_aa: return True
    else: return False
    
    
    
def missing_aa_inside_protein(biopython_structure):
    '''check if any missing/unsolved residues in structure in the 
    middle of the sequence. Checks based on structure not on jumps
    in sequence (the latter can happen when lab cuts out part of seq)
    ''' 
    
    ppb = PPBuilder()
    for i, chain in enumerate(biopython_structure):
        pps = ppb.build_peptides(chain)
        if len(pps) > 1:
            # chain breaks in structure (not in sequence) if missing residues in pps => pps>1
            return True

        elif len(pps)==1:
            return False


    
def too_big(mmtf_record, max_len):
    size = get_size(mmtf_record)
    if size > max_len: return True
    else: return False
    
    
def too_small(mmtf_record, min_len):
    size = get_size(mmtf_record)
    if size < min_len: return True
    else: return False


def get_size(mmtf_record):
    '''get total residue size of all chains in protein''' 

    protein_entities = [entity for entity in  mmtf_record.entity_list if entity['type'] == 'polymer']
    size = 0
    for chain in protein_entities:
        sequence = chain['sequence']
#         if len(sequence) < 50:
#             print(f"OBS! protein has chain with nr res<50 - reconsider.\n {mmtf_record.title}, len:{len(sequence)}")
        size+=len(sequence)
    
    return size

def get_biopython_structure(mmtf_record):
    structure = BioMMTF.get_from_decoded(mmtf_record)
    first_model = structure.get_list()[0]
    return first_model


def unsuitable_ligands(mmtf_record):

    # define if ligand is unsuitable/not allowed
    unsuitable = False

    # get all ligands/non_proteins in pdb
    non_proteins = [entity for entity in  mmtf_record.entity_list if entity['type'] == 'non-polymer']

    # if len(non_proteins) > 0: # if = 0 => no other molecyles than protein
    for entity in non_proteins:
        ligand = entity['description']
        if ligand not in allowed_metals and ligand not in standard_crystallization_ligands:
            unsuitable = True

    return unsuitable  
                
                

              




