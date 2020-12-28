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
traj_format = 'xtc'

#with open('CATH4.2_cleaned_data_no_NMR_res2.0_noMetals.dat','r') as f:
#    cleaned_pdbs = f.readlines()
with open('NMRs_rmsd_less3.5.txt','r') as f:
    cleaned_pdbs = f.readlines()
cleaned_pdbs = [pdb.strip() for pdb in cleaned_pdbs]
print(len(cleaned_pdbs))
# filter out nmr anyway

count_already_sim = 0

# test run 5 proteins
#for name_pdb in cleaned_pdbs[200:500]: # bio assemblies
for name_pdb in cleaned_pdbs[:85]: # NMRs
    #path_pdb = f'{cwd}/data/pdbs_cleaned_bio_assembly/{name_pdb}_clean.pdb'
    #out_dir = f'{cwd}/simulations_new/{name_pdb}/'
    path_pdb = f'{cwd}/data/pdbs_cleaned_NMRs/{name_pdb}_clean.pdb'
    out_dir = f'{cwd}/simulations_new/NMRs/{name_pdb}/'
    cmd = f'{MD_python} {cwd}/protocol_prototype_read_from_pdbs.py --pdb {path_pdb} --trajectory_format {traj_format} -o_dir {out_dir} --cuda'
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    with open(out_dir + submit_file_name, 'w') as submit:
        for line in submit_file_cont:
            submit.write(line+'\n')
        submit.write(f'#PBS -N NMR_{name_pdb}\n')
        submit.write(f'#PBS -o log_sim_{name_pdb}.log\n')
        submit.write('\n')
        submit.write(cmd)
    
    
    print(name_pdb, 'finish')
    
    # submit from seperate pdb dir
    os.chdir(out_dir)
    os.system(f'qsub {submit_file_name}')
    os.chdir(cwd)

    
    
    
    
    
