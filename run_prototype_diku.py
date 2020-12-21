import glob, os



submit_file_name = 'submit_diku.sh'
submit_file_cont = ['#!/bin/bash', 
                    '#SBATCH --ntasks=1',
                    '#SBATCH --cpus-per-task=5',
                    '#SBATCH --time=05:00:00',
                    '#SBATCH -p boomsma',
                    '#SBATCH --gres=gpu:1']

cwd = os.getcwd() 
MD_python = '/home/trz846/anaconda3/envs/MD/bin/python'
traj_format = 'xtc'

# remove already simulated
with open('list_cleaned_pdbs.txt','r') as f:
    cleaned_pdbs = f.readlines()
cleaned_pdbs = [pdb.strip() for pdb in cleaned_pdbs]
with open('already_sim_diku_comp.txt','r') as f:
    already = f.readlines()
already = [pdb.strip() for pdb in already]
cleaned_pdbs = [pdb for pdb in cleaned_pdbs if pdb not in already]


# test run 5 proteins
for name_pdb in cleaned_pdbs[400:402]:
    path_pdb = f'{cwd}/data/pdbs_cleaned/{name_pdb}_clean.pdb'
    out_dir = f'{cwd}/simulations/{name_pdb}/'
    cmd = f'{MD_python} {cwd}/protocol_prototype.py --pdb {cwd}/{path_pdb} --trajectory_format {traj_format} -o_dir {out_dir} --cuda'
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

    
    
    
    
    
