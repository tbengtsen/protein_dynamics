'''
Cleans raw downloaded pdbs to parse without errors/crashing for protocol. 
'''
import sys, glob, os 
sys.path.insert(0, '../scripts/')
from clean_pdb import *

# reduce executable
reduce_executable = '../scripts/reduce/reduce'
reduce_het_dict = '../scripts/reduce/reduce_wwPDB_het_dict.txt'
out_dir = './pdbs_cleaned/'
raw_pdbs = glob.glob('pdbs_raw/*.pdb')
log_file = 'log_clean_raw_pdbs.log'

assert len(raw_pdbs)!= 0, ('Cannot find any downloaded pdbs.',
                           ' Please run: download_raw_pdbs.py first ')

assert os.path.isfile(reduce_executable), (
        f'Cannot find reduce executable in scripts',
        'please download. See installation guide in README.md')


assert os.path.isfile(reduce_het_dict), (
        f'Cannot find reduce file reduce_wwPDB_het_dict.txt',
        'please download and put in same directory as reduce executable')


# settings to run 
log = open(log_file, 'w')


for pdb_file in raw_pdbs:
    try:
        clean_pdb(pdb_file, out_dir, reduce_executable)
    except:
        log.write(f'{pdb_file}, something went wrong when cleaning the file')
        
    
log.close()