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
