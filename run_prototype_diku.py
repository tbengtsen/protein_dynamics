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
path_cleaned = 'data/CATH4.3/pdbs_cleaned/' ## XX <---- OBS ! CHANGE
cleaned_pdbs = glob.glob(path_cleaned +'/*.pdb')
traj_format = 'xtc'

# test run 5 proteins
for path_pdb in cleaned_pdbs[50:75]:
    name_pdb = path_pdb.split('/')[-1].split('_')[0]
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

    
    
    
    
    
