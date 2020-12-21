import glob, os

submit_file_name = 'submit_computerome.sh'
submit_file_cont = ['#!/bin/bash',
        '#PBS -l nodes=1:ppn=5:gpus=1',
                    '#PBS -l mem=24gb',
                    '#PBS -l walltime=12:00:00',
                    '#PBS -q batch',
                    '#PBS -W group_list=ku_00043',
                    '#PBS -A ku_00043']




# paths 
MD_python = '/home/projects/ku_00043/people/toneben/miniconda3/envs/MD/bin/python'
cwd = os.getcwd()
path_cleaned = 'data/pdbs_cleaned/'  
cleaned_pdbs = glob.glob(path_cleaned +'/*.pdb')
traj_format = 'xtc'

# test run 5 proteins
for path_pdb in cleaned_pdbs[575:900]:
    name_pdb = path_pdb.split('/')[-1].split('_')[0]
    out_dir = f'{cwd}/simulations/{name_pdb}/'
    cmd = f'{MD_python} {cwd}/protocol_prototype.py --pdb {cwd}/{path_pdb} --trajectory_format {traj_format} -o_dir {out_dir} --cuda'
    print('\ncmd:\n',cmd)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    with open(out_dir + submit_file_name, 'w') as submit:
        for line in submit_file_cont:
            submit.write(line+'\n')
        submit.write(f'#PBS -N {name_pdb}\n')    
        submit.write(f'#PBS -o log_sim_{name_pdb}.log')
        submit.write('\n')
        submit.write(cmd)
    
    print(name_pdb, 'finish')
    
    # submit from seperate pdb dir
    os.chdir(out_dir)
    os.system(f'qsub {submit_file_name}')
    os.chdir(cwd)

    
    
    
    
    
