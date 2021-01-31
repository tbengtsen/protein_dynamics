"""
Class that sets up object to contain protein g public data. 
For ease of use when benchmarking stability predictors

"""
import json, os, sys, glob, time, copy
import numpy as np 
import pandas as pd


class Allstability:
    '''make stability object '''
    def __init__(self, data_type = None):

        self.main_path = '../Data/'
        
        if not os.path.exists(self.main_path): sys.exit('Cannot find main_path to Data directory')
        self.protein_features = self.read_json(self.main_path  + 'protein_features.json')
        self.protein_info = self.read_json(self.main_path  + 'info_all_stability.json')
        self.ddG = self.read_stability(self.main_path + 'Processed/ddG/')
        self.proteins = self.get_all_proteins()
        print(f'\n\n\n**Contains ddG single mutations for the following proteins:**\n\n{self.ddG.keys()}')
        print(f"\n\nRead about the protein data and the associated readme files for each protein by: object.protein_info[name] or object[protein_name].README\n")
      
    
    
    def __getitem__(self, name):
        '''returns dictionay with all data associated with the given protein name'''
        
        
        data = {}
        
        # check exist 
        if name not in self.protein_info.keys():
            sys.exit(f'{name} not in stability data class')
        
        # add basic info
        data['protein_info'] = self.protein_info[name]
        data['README'] = self.protein_info[name]['README']
        
        
        # add ddG   
        if name in self.ddG.keys():   
            data['ddG'] = self.ddG[name]
        else:
            data['ddG'] = None

## ==== load Mayo's protein G data

 import json, os, sys, glob, copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

sys.path.insert(0, '/home/trz846/stability_pipeline/Scripts/')
from all_stability import AllStability

#Initialise/load stability object with stability data from Mayos protein G public data
stability = AllStability()

# c = 'red'
# ## SURFACE CORRELATION TO ddGS ##
nr_proteins = len(stability.ddG.keys())
fig, ax = plt.subplots(1,nr_proteins,sharey=False, constrained_layout=True,figsize=(25,4))
fig.suptitle('Experiment -  min ddG per position vs Secondary structure \n  increased stability <-----',fontsize=25)

plot = [0]

for protein in stability.ddG.keys():
    try:
        sec_struc = stability[protein]['Secondary_structure']
        ddG = stability.ddG[protein]
        sec_struc_column = [sec_struc [pos-1] for pos in ddG['position']] #index - position-1
        ddG['sec_struc'] = sec_struc_column
        ddG.sort_values(by='sec_struc')
        ddG_corr = ddG.dropna() # remove any NaN
        # spearman  = stats.spearmanr(ddG_corr['ddG'], ddG_corr['SASA'])[0] direct correlation

        # plot correlation betweeen SASA and maximum improvement factor in each position in stability
        min_ddG = ddG_corr.groupby(['position']).min() # max of each position



# =========================
# visualise in NGLVIEW
# =========================
protein_name = 'protein_g'

# --- ball and stick visualisation ---
import nglview as nv
import MDAnalysis
exp = 'ddg'
df = stability.ddg[protein_name] # <----- <-----INPUT <------ <------  ¯\_(ツ)


top_hits_exp = df.sort_values(exp) [:10]
resid_exp = top_hits_exp['position']
resid_exp = ' '.join([str(res) for res in resid_exp])



pdb_file = glob.glob(f'../Data/PDBs/{protein_name}/*.pdb')[0]

# load pdb into mdanalysis
u = MDAnalysis.Universe(pdb_file ,pdb_file )


w = nv.show_mdanalysis(u,gui=True)
w.clear_representations(component=0)
w.add_representation('cartoon', color='blue',opacity=0.5, component=0)
w.add_representation("ball+stick",selection=resid_exp,color='red', component=0)

w

# ---- Color after hits ---
import math
protein_name = 'protein_g' # <----- <-----INPUT <------ <------  ¯\_(ツ)
exp = 'ddg' 
ddg = stability.ddg[protein_name] # <----- <-----INPUT <------ <------  ¯\_(ツ)
cutoff = 1 # <----- <-----INPUT PLAY AROUND WITH THIS TO MAKE CONTINUES COLORSCALE <------ <------  ¯\_(ツ)


max_ddg = ddg.groupby(['position']).max() # max of each position


pdb_file = glob.glob(f'../Data/PDBs/{protein_name}/*.pdb')[0]

# load pdb into mdanalysis
u = MDAnalysis.Universe(pdb_file ,pdb_file )
u.add_TopologyAttr('tempfactors') # add empty attribute for all atoms
protein = u.select_atoms('protein and name CA') # select protein atoms

# add values where missing 
ddg_present = []
ddg_coloring = []
ddg_missing = []
for residue in protein.residues:
    try:
        if max_ddg['ddg'][residue.resid] > cutoff:
            ddg_coloring.append(math.log(max_ddg['ddg'][residue.resid]))
        elif max_ddg['ddg'][residue.resid] < cutoff:
            ddg_coloring.append(math.log(cutoff))
        ddg_present.append(str(residue.resid))
    except:
        ddg_missing.append(str(residue.resid))
        ddg_coloring.append(math.log(cutoff))

    
# add b factor/tempfactor to residues
for residue, ddg_value in zip(protein.residues, ddg_coloring):
        residue.atoms.tempfactors = -ddg_value 
# find top hits to visualise as ball+stick
top_hits = ddg.sort_values(by='ddg',ascending=False)
top_hits =  top_hits ['position'][:5]
top_hits = [str(pos) for pos in top_hits]
top_hits = ' '.join(top_hits)

w = nv.show_mdanalysis(u,gui=True)
w.clear_representations(component=0)
w.add_representation('cartoon',color='cyan', component=0,opacity=0.3)
w.add_representation("ball+stick", selection=top_hits, color='blue', component=0)
w.add_representation('cartoon', selection=' '.join(ddg_missing), color='black', component=0)
w.add_representation('backbone', selection=top_hits,color ='blue', component=0)

w.add_representation('cartoon', selection=' '.join(ddg_present), color='bfactor', component=0)

w


