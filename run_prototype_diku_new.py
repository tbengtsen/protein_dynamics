import glob, os



submit_file_name = 'submit_diku.sh'
submit_file_cont = ['#!/bin/bash', 
                    '#SBATCH --ntasks=1',
                    '#SBATCH --cpus-per-task=5',
                    '#SBATCH --time=04:00:00',
                    '#SBATCH -p gpu',
                    '#SBATCH --gres=gpu:1']

cwd = os.getcwd() 
MD_python = '/home/trz846/anaconda3/envs/MD/bin/python'
traj_format = 'xtc'

# remove already simulated
with open('CATH4.2_cleaned_data_no_NMR_res2.0_noMetals.dat','r') as f:
    cleaned_pdbs = f.readlines()
cleaned_pdbs = [pdb.strip() for pdb in cleaned_pdbs]
print(len(cleaned_pdbs))
# filter out nmr anyway
with open('NMR_pdb_ids.txt','r') as f:
    NMRs = f.readlines()
NMRs = [nmr.strip() for nmr in NMRs]

count_nmr = 0
count_already_sim = 0

# test run 5 proteins
for name_pdb in cleaned_pdbs[:250]:
    if os.path.isdir(f'{cwd}/simulations_new/{name_pdb}'):
        print(f'{name_pdb} already simulated')
        count_already_sim +=1
        continue
    else:    
        path_pdb = f'{cwd}/data/pdbs_cleaned_bio_assembly/{name_pdb}_clean.pdb'
        out_dir = f'{cwd}/simulations_new/{name_pdb}/'
        cmd = f'{MD_python} {cwd}/protocol_prototype_read_from_pdbs.py --pdb {path_pdb} --trajectory_format {traj_format} -o_dir {out_dir} --cuda'
        print('\ncmd:\n',cmd)
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        with open(out_dir + submit_file_name, 'w') as submit:
            for line in submit_file_cont:
                submit.write(line + '\n')
            submit.write(f'#SBATCH --job-name={name_pdb}\n')
            submit.write('\n')
            submit.write(cmd)

        print(name_pdb, 'finish')

        # submit from seperate pdb dir
        os.chdir(out_dir)
        os.system(f'sbatch {submit_file_name}')
        os.chdir(cwd)

print(f'discarded {count_nmr} pdbs due to nmr')
print(f'discarded {count_already_sim} due to already sim')
    
    
    
    
