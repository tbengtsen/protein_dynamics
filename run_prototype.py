import glob, os, sys, shutil



submit_file = 'submit_computerome.sh'
MD_python = '/home/projects/ku_00043/people/toneben/miniconda3/envs/MD/bin/python'

path = 'data/pdbs_cleaned/'
cleaned_pdbs = glob.glob(path+'/*.pdb')
traj_format = 'xtc'

# test run 5 proteins
for path_pdb in cleaned_pdbs[:2]:
    name_pdb = path_pdb.split('/')[-1].split('_')[0]
    out_dir = 'simulations/' + name_pdb + '/'
    cmd = f'{MD_python} run_prototype.py --pdb {path_pdb} --trajectory_format {traj_format} --out_dir {out_dir}'
    print('\ncmd:\n',cmd)
    if not os.path.isdir(out_dir):
        os.makedirs('out_dir')
    os.popen(f'cp {submit_file} {out_dir}')
    with open(out_dir + submit_file, 'a+') as submit:
        submit.write('\n')
        submit.write(cmd)
    
    print(name_pdb, 'finish')
    

    
    
    
    
    
