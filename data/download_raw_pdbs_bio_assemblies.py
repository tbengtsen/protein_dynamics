import sys, os,json
import urllib.request      


# load all PDB ids from CATH4.3 that are possible to simulate
#clean_data_path = 'CATH4.2/CATH4.2_cleaned_data.json'
#clean_data_path = 'CATH4.2/CATH4.2_cleaned_data_no_NMR_res2.5.json'
#with open(clean_data_path,'r') as f:
#    clean_data = json.load(f)
with open('../CATH4.2_cleaned_data_no_NMR_res2.0_noMetals.dat', 'r') as f:
    clean_data = f.readlines()
clean_data = [pdb.strip() for pdb in clean_data]


# define download directory
target_directory = 'pdbs_raw_bio_assembly/'
target_path = os.path.dirname(target_directory)
if not os.path.exists(target_path):
    os.makedirs(target_path)
count = 0
# download all pdbs from processed CATH4.3 from rcsb 
for pdb  in clean_data:
    target_file = target_path + f'/{pdb}.pdb'
    if not os.path.isfile(target_file):
        response = None
        try:
            os.system(f'wget ftp://ftp.wwpdb.org/pub/pdb/data/biounit/PDB/all/{pdb}.pdb1.gz {target_path}')
        except:
            pass
        try:
            os.system(f'wget ftp://ftp.wwpdb.org/pub/pdb/data/biounit/PDB/all/{pdb}.pdb2.gz {target_path}')
        except:
            pass
        try:
            os.system(f'wget ftp://ftp.wwpdb.org/pub/pdb/data/biounit/PDB/all/{pdb}.pdb3.gz {target_path}')
        except:
            pass

