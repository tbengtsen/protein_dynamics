'''
Cleans raw downloaded pdbs to parse without errors/crashing. 
'''
import sys, glob
sys.path.insert(0, '../scripts/')
from clean_pdb import *


# settings to run 

raw_pdbs = glob.glob('pdbs_raw/*.pdb')

for pdb_file in raw_pdbs:
    
