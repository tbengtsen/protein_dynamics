import glob, os



submit_file_name = 'submit_diku.sh'
submit_file_cont = ['#!/bin/bash', 
                    '#SBATCH --ntasks=1',
                    '#SBATCH --cpus-per-task=5',
                    '#SBATCH --time=8:00:00',
                    '#SBATCH -p gpu',
                    '#SBATCH --exclude=a00863,a00862,a00861,a00860',
                    '#SBATCH --gres=gpu:1']

cwd = os.getcwd() 
MD_python = '/home/trz846/anaconda3/envs/MD/bin/python'
traj_format = 'xtc'

with open('empty_run_files.dat','r') as f:
    cleaned_pdbs = f.readlines()
cleaned_pdbs = [pdb.strip() for pdb in cleaned_pdbs]
print(len(cleaned_pdbs))

for name_pdb in cleaned_pdbs: 
    if os.path.isdir(f'{cwd}/simulations_new/{name_pdb}'):
        print(f'{name_pdb} already simulated')
        continue
    else:    
        
        path_pdb = f'{cwd}/data/pdbs_cleaned_new/{name_pdb}_clean.pdb'
        out_dir = f'{cwd}/simulations_new/{name_pdb}/'
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        cmd = f'{MD_python} {cwd}/protocol_prototype_read_from_pdbs.py --pdb {path_pdb} --trajectory_format {traj_format} -o_dir {out_dir} --cuda'
        
        with open(out_dir + submit_file_name, 'w') as submit:
            for line in submit_file_cont:
                submit.write(line + '\n')
            submit.write(f'#SBATCH --job-name=er_{name_pdb}\n')
            submit.write('\n')
            submit.write(cmd)

        print(name_pdb, 'finish')

        # submit from seperate pdb dir
        os.chdir(out_dir)
        os.system(f'sbatch {submit_file_name}')
        os.chdir(cwd)

    
    
    
    
