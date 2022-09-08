'''

This file splits the all simulated pdb data into train, val, and test and saves it in a json
for use in J. Ingrahams Graph transformer [1]. 

The datasplit itself are taken from Ingrahams datasplit to be comparable, see citation below [1].

The notebook allows for two types of split:

   1) splits with simulations

   2) splits with single pdbs (original pdbs) typical for reference, see flag options
   
[1] author = {Ingraham, John and Garg, Vikas K and Barzilay, Regina and Jaakkola, Tommi},
title = {Generative Models for Graph-Based Protein Design},
booktitle = {Advances in Neural Information Processing Systems}
year = {2019}


'''

import json, random, os, argparse, warnings


def parser():
    '''
    parser for script
    '''
    
    parser = argparse.ArgumentParser(description='Creates data split container in json format for use with J. Ingrahams Graph Transformer')
    parser.add_argument('-list_pdbs',
                        dest='list_pdbs',
                        type=str,  
                        required = False,
                        default='../data/CATH4.2_ALL_STATUS.json',
                        help='list or json that contains the pdbs names of simulations split by. Default ../data/CATH4.2_ALL_STATUS.json ')
    parser.add_argument('-ingraham_split',
                        dest='ingr_split',
                        type=str,  
                        required = False,
                        default = '../data/Ingrahams_chain_set_splits.json',
                        help='path to json file from Ingraham with data split. Default: ../data/Ingrahams_chain_set_splits.json')
    
    parser.add_argument('-o',
                        dest='out', 
                        type=str,
                        required = False,
                        default='./data_split_{}ps.json',
                        help='output json file with data split. Default: ./data_split_{freq}ps.json')
    
    parser.add_argument('-freq',
                        dest='freq',
                        type=int,
                        default = 1000,
                        required = False,
                        help='how often to use output frames from simulation in ps. Default = 1000 ps' )
    
    parser.add_argument('-incl_orig_pdb',
                        dest='orig_pdb',
                        action="store_true",
                        help='include structure of original pdb from PDB database into split data in addition to simulation frames. Default False' )
    
    parser.add_argument('-max_len',
                        dest='max_len',
                        type=int,
                        default = 500,
                        required = False,
                        help='Filter function. Only use simulations of proteins <=max_len in all splits. Default <=500' )
    
    parser.add_argument('-exp_method',
                        dest='exp',
                        type=str,
                        default = 'all',
                        choices = ['Xray','NMR','all'],
                        required = False,
                        help='Filter function. Only use simulations of proteins experimental solved by either Xray crystallization, NMR, or all Default:all' )
    
    parser.add_argument('-max_chains',
                        dest='max_chains',
                        type=int,
                        default = -1,
                        required = False,
                        help='Filter function. Set the maximum number of chains allowed in the protein. Default: all chains (-1)' )
    

    parser.add_argument('-max_pdbs_in_train',
                        dest='max_train',
                        type=int,
                        default= -1 ,
                        required = False,
                        help='Filter function. Only use a specific number of pdbs in the train dataset. Default: All  ' )

    parser.add_argument('-output_reference_single_split',
                        dest='out_singles',    
                        action="store_true",
                        help='whether to output a split file for single structures to use for reference.' )


    
    args = parser.parse_args()

    
    return args




def read_ingraham_split(path):
    '''Reads the original Ingraham data split file to use for splitting
    '''
    
    with open(path,'r') as f:
        ingr_split = json.load(f)
    train_ingr = ingr_split['train']
    test_ingr = ingr_split['test']
    val_ingr = ingr_split['validation']
    
    
    return train_ingr, test_ingr, val_ingr




def read_sim_status(path,
                    filter_max_len = None,
                    filter_exp_method = None,
                    filter_max_chains = -1,
                    filter_max_train = -1,):
    '''
    Reads the CATH4.2_ALL_STATUS.json and returns pdbs which have finished
    their simulations. 
    Additionally filters out either by lenght of proteins, by experimental 
    method or by nr of chains.
    '''
    

    with open(path,'r') as f:
            all_cath = json.load(f)
    
    
    # get all pdbs that have been simulated
    simulated_pdbs = [pdb for pdb in all_cath if all_cath[pdb]['simulation status'] == 'Finished' ]

    # filter out pdbs based on filter conditions
    if  filter_max_len != None:
        simulated_pdbs  = filter_pdbs_by_len(all_cath, simulated_pdbs, filter_max_len)
        
    if filter_exp_method != None and filter_exp_method != 'all' :
        simulated_pdbs  = filter_pdbs_by_exp(all_cath, simulated_pdbs, filter_exp_method)
        
    if filter_max_chains != -1 :
        simulated_pdbs  = filter_pdbs_by_chains(all_cath, simulated_pdbs, filter_max_chains)
        
    
    # split into datasplit based on Ingraham's datasplit
    train_sim = [pdb for pdb in simulated_pdbs if all_cath[pdb]['split']=='train']
    test_sim = [pdb for pdb in simulated_pdbs if  all_cath[pdb]['split']=='test']
    val_sim = [pdb for pdb in simulated_pdbs if all_cath[pdb]['split']=='validation']
    
    # filter by max number of train pdbs 
    if filter_max_train != -1:
         train_sim = filter_pdbs_by_max_train(all_cath, train_sim, filter_max_train)
    
    
    print(f'\nSize of train datasplit: {len(train_sim)}')
    print(f'Size of test datasplit: {len(test_sim)}')
    print(f'Size of val datasplit: {len(val_sim)}\n')
    
    # ensure that no datasplit are empty given filter conditions
    for datasplit in [train_sim,test_sim,val_sim]:
        assert len(datasplit) > 0, f' The choosen set of filterings leaves no proteins in one of the datasplits'
    
    

    return train_sim, test_sim, val_sim




def filter_pdbs_by_len(all_cath, simulated_pdbs, filter_max_len):
    '''Filter out all proteins that have more residues than filter_max_len'''
    
    simulated_pdbs = [pdb for pdb in simulated_pdbs if all_cath[pdb]['size'] <= filter_max_len]
    
    return simulated_pdbs 


def filter_pdbs_by_exp(all_cath, simulated_pdbs, filter_exp_method):
    '''filter (keep only) list of pdbs by experimental method'''
    
    if filter_exp_method == 'Xray': method = 'X-RAY DIFFRACTION'  
    if filter_exp_method == 'NMR' : method = 'SOLUTION NMR' 
    simulated_pdbs = [pdb for pdb in simulated_pdbs if method in all_cath[pdb]['method']]

    return simulated_pdbs


def filter_pdbs_by_chains(all_cath, simulated_pdbs, filter_max_chains):
    '''Filter out all proteins with nr chains > filter_max_chains'''
    
    simulated_pdbs = [pdb for pdb in simulated_pdbs if all_cath[pdb]['chains'] <= filter_max_chains]
    
    return simulated_pdbs 

def filter_pdbs_by_max_train(all_cath, simulated_pdbs, filter_max_train):
    '''remove excess pdbs in training dataset if a certain number of pdbs
    are specified for training. 
    '''
    
    simulated_pdbs = simulated_pdbs[:filter_max_train]
    
    if filter_max_train > len(simulated_pdbs):
        warnings.warn('The maximum number of pdbs specified in train is higher than number of the number of finished simulations ')
    
    return simulated_pdbs



def get_data_split (split_ingraham:list, 
                    sim:list,
                    freq:int, 
                    incl_orig_pdb = True):

    '''
    Split the pdbs which have been simulated after the Ingraham data split. 
    Returns dictionary with splits both for frames in simulations and one 
    for the corresponding single pdb structures to compare with.
    '''

    split_simulations = []
    split_orig_pdbs = []
    
    for pdb_chain in split_ingraham : # uses split ingraham to get chains as well as pdb id
        pdb = pdb_chain.split('.')[0]
        
        if pdb in sim:
            
            # add pdb to single_reference like so: pdb.chain 
            split_orig_pdbs.append(pdb_chain)
            
            ## add sim frames to pdb+chain like so: pdb.chain.frame
            
            # append original pdb denoted -1 
            if incl_orig_pdb:
                pdb_frame = pdb_chain + '.-1'
                split_simulations.append(pdb_frame)
            
            # append first frame in sim = 50 ps 
            # as did not save sim starting structure (the one after equil)
            if freq !=50:
                pdb_frame = pdb_chain + f'.50'
                split_simulations.append(pdb_frame) 
            
            # append all rest frames in structure with freq intervals
            for frame in range(freq, 20000, freq): # pick out frames with freq up to 20ns
                pdb_frame = pdb_chain + f'.{frame}'
                split_simulations.append(pdb_frame)

            # append last frame (19950ps) - bc frame 20000ps not saved in sim
            pdb_frame = pdb_chain + f'.19950'
            split_simulations.append(pdb_frame)
                

    return split_simulations, split_orig_pdbs




if __name__=='__main__':
    
    args = parser()
    
    # check input arguments 
    assert os.path.exists(args.list_pdbs), 'the input list or json given to -list_pdbs does not exists '
    assert os.path.exists(args.ingr_split), 'the path to the ingraham data split file does not exists'

    
    train_sim, test_sim, val_sim = read_sim_status(
                                        args.list_pdbs,
                                        filter_max_len = args.max_len,
                                        filter_exp_method = args.exp,
                                        filter_max_chains = args.max_chains,
                                        filter_max_train  = args.max_train)
    
    train_ingr, test_ingr, val_ingr = read_ingraham_split(args.ingr_split)
    
    # keep same test and validation test for all trained number of pdbs:
    sim_split = {'train':[],'test':[],'validation':[]}
    single_split  = {'train':[],'test':[],'validation':[]} 
    
    ## gather all splits with simulations 
    sim_split['test'], single_split['test'] = get_data_split (test_ingr, 
                    test_sim,
                    freq = args.freq,
                    incl_orig_pdb = args.orig_pdb)
    
    sim_split['validation'], single_split['validation'] = get_data_split (val_ingr, 
                    val_sim,
                    freq = args.freq,
                    incl_orig_pdb = args.orig_pdb)
    
    # use all simulated train pdbs
    sim_split['train'], single_split['train'] = get_data_split (train_ingr, 
            train_sim,
            freq = args.freq,
            incl_orig_pdb = args.orig_pdb)
    
        
    # save in json
    try:
        o_name = f'{args.out}'.format(args.freq) 
    except:
        o_name = args.out
    
    # add file type if missing
    if args.out.split('.')[-1] == 'json':
        o_name = args.out
    else:
        o_name = args.out + '.json'
    
    # save simulations split 
    with open(o_name, 'w') as f:
        json.dump(sim_split,f)
    
    # save single split (here chains in included to use with json chain_split.jsonl)
    if args.out_singles:
        print(o_name)
        # also output the same splits with single original data referencing
        o_name_singles = o_name.split('.json')[0] + '_single_structures.json'
        with open(o_name_singles, 'w') as f:
            json.dump(single_split,f)    


   