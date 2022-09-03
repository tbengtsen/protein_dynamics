import sys, os,json
import urllib.request      


# load all PDB ids from CATH4.2 that are appropiate for simulations
clean_data_path = 'CATH4.2/CATH4.2_cleaned_data.json'
with open(clean_data_path,'r') as f:
    clean_data = json.load(f)


# define download directory
target_directory = 'pdbs_raw/'
target_path = os.path.dirname(target_directory)
if not os.path.exists(target_path):
    os.makedirs(target_path)

# download all pdbs from processed CATH4.2 from rcsb 
for ds_type  in clean_data:
    for pdb in clean_data[ds_type]:
        target_file = target_path + f'/{pdb}.pdb'
        if not os.path.isfile(target_file):
            try:
                response = urllib.request.urlopen(f'https://files.rcsb.org/download/{pdb}.pdb')
                with open(target_file, 'wb') as f:
                    f.write(response.read())
                print(f'downloaded {pdb}')    

            except: 
                print(f'cannot download {pdb}')

