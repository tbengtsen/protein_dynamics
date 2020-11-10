import sys, os,json
import urllib.request      


# # load all PDB ids from CATH4.3 that are possible to simulate
# clean_data_path = 'CATH4.3/CATH_cleaned_data.json'
# with open(clean_data_path,'r') as f:
#     clean_data = json.load(f)


clean_data = ['2lzw']
# define download directory
target_directory = 'pdbs_raw/'
target_path = os.path.dirname(target_directory)
if not os.path.exists(target_path):
    os.makedirs(target_path)




for pdb in clean_data:
    try:
        response = urllib.request.urlopen(f'https://files.rcsb.org/download/{pdb}.pdb')
        target_file = target_path + f'/{pdb}.pdb'

        with open(target_file, 'wb') as f:
            f.write(response.read())
            
    except: 
        print(f'cannot download {pdb}')
        