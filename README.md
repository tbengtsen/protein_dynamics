# protein_dynamics
Dynamics augmented protein structure dataset 




# Installation
Install reduce. This program is used by the parser to set charges and add missing hydrogens to the proteins
```
cd scripts/
git clone https://github.com/rlabduke/reduce.git 
cd reduce/
make; make install # This might give an error but provide the reduce executable in this directory
```



# Data/simulation tracker
The data used for the simulations performed here are based on the CATH4.2
dataset[1], see https://www.cathdb.info/wiki.  

```
# how to retrive data on simulation tracker
with open('data/CATH4.2_ALL_STATUS.json', 'r') as f:
    all_cath42 = json.load(f)
    

# get all pdbs in CATH4.2
pdbs = list(all_cath42.keys())

# print a single pdb to see structure of tracker dictionary
print(all_cath42['3k4u'])

  {'split': 'train',
 'prefilter status': 'removed',
 'prefilter reasons': ['low resolution', 'gap in chain', 'unsuitable ligand',  'too big'],
 'cath_nodes': ['3.40.190'],
 'size': 1521,
 'chains': 6,
 'method': ['X-RAY DIFFRACTION'],
 'simulation status': 'N.A.',
 'simulation reason': 'N.A.',
 'rmsd NMR': nan,
 'resolution': 2.62}

# too see if simulation of protein is performed: 
simulated_pdbs = [pdb for pdb in all_cath42 if all_cath42[pdb]['simulation status'] == 'Finished']
```


# Visualise a single simulation
see: 
```scripts/notebook_protocol/view_single_simulation.ipynb```



# Citations:
```
1) 
@article{CATH,
author = {Sillitoe, Ian and Bordin, Nicola and Dawson, Natalie and Waman, Vaishali P and Ashford, Paul and Scholes, Harry M and Pang, Camilla S M and Woodridge, Laurel and Rauer, Clemens and Sen, Neeladri and Abbasian, Mahnaz and Le Cornu, Sean and Lam, Su Datt and Berka, Karel and Varekova, Ivana Hutařová and Svobodova, Radka and Lees, Jon and Orengo, Christine A},
title = "{CATH: increased structural coverage of functional space}",
journal = {Nucleic Acids Research},
volume = {49},
number = {D1},
pages = {D266-D273},
year = {2020},
doi = {10.1093/nar/gkaa1079}}
```

```
2) 
@inproceedings{ingraham2019generative,
author = {Ingraham, John and Garg, Vikas K and Barzilay, Regina and Jaakkola, Tommi},
title = {Generative Models for Graph-Based Protein Design},
booktitle = {Advances in Neural Information Processing Systems}
year = {2019}
}
```