This directory contains files used with the CATH4.2 protein structure classification database [1].
to filter out PDBs not appropiate for this type of automated simulations
# Overview of files and their function

*cleaned_data.json*
Overview of PDB filtering. Contains all protein PDB IDs found CATH4.2, their individual basic structural info, and whether they are deemed appropiate for basic automated simulations based on filtering criteria lined out in the file 'filter_cath_dataset.ipynb', as well as their individual data set split, i.e. train, test, or validation. The data split are copied from [2] for later comparison to results
Example: 

'''
{"2fyz": {"split": "train", "prefilter reasons": ["too big"], "cath_nodes": ["1.20.5"], "size": 539, "chains": 6, "method": ["X-RAY DIFFRACTION"], "prefilter status": "removed", "simulation status": "N.A.", "simulation reason": "N.A.", "rmsd NMR": NaN, "resolution": 2.2}, "4nav": {"split": "train", "prefilter reasons": ["low resolution", "too big"], "cath_nodes": ["3.40.50"], "size": 763, "chains": 4, "method": ["X-RAY DIFFRACTION"], "prefilter status": "removed", "simulation status": "N.A.", "simulation reason": "N.A.", "rmsd NMR": NaN, "resolution": 2.69},
'''

*Ingrahams_chain_set_splits.json*
Data split of all chains from CATH4.2 PDB IDs as used in [2] 

*filter_cath_dataset.ipynb*
filtering file which loops through all PDB IDs in CATH4.2 and obtainsor calculates structural information about the PDB to filter out PDBs deemed inappropiate for basic automated simulations based on filtering criterias lined out in the file

*membrane_protein_database_mpstruc.xml*
File used to filter out PDB IDs that are membrane PDBs and thus inappropiate for this type of automated simulations


# Citations:
```
1)
@article{CATH,
author = {Sillitoe, Ian and Bordin, Nicola and Dawson, Natalie and Waman, Vaishali P and Ashford, Paul and Scholes,
Harry M and Pang, Camilla S M and Woodridge, Laurel and Rauer, Clemens and Sen, Neeladri and Abbasian, Mahnaz and Le
<C2><A0>Cornu, Sean and Lam, Su Datt and Berka, Karel and Varekova, Ivana<C2><A0>Huta<C5><99>ov<C3><A1> and Svobodov
a, Radka and Lees, Jon and Orengo, Christine A},
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
}`
