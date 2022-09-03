# DATA INFORMATION
The PDB data used for these simulations are based on the CATH4.2 [1] Protein Structure Classification database which filters and classifies protein domains based on evolutionary relationship. This is used to obtain a balanced dataset. The dataset split used are copied from John Ingraham's paper [2] for later comparison of results. However, here we filter out more proteins based on whether their structure are appropiate for these types of automated simulations. 

# Simulation filter criterias
The conditions used to filter out proteins are as follows: 

    - is a membrane protein
    - unsuitable experimental method  for solving the structure (unsuitable only in the sense of usefull for simulations)
    - low X-ray resolution 
    - to high RMSD (>3.5) between models in NMR models 
    - contains unknown amino acids
    - gaps/unresolved parts inside the chain of a protein (in terminals these are allowed)
    - unsuitable ligand(s) (unsuitable only in the sense of usefull for automation of simulations)
    - too big a protein (>500 residues)
    - any chain too small in the  protein (<50 residues). Meaning it is unlikely to process a (somewhat) 3D structure by itself

**OBS!** It should be stressed that neither of the filtering conditions descriped above necessarily say anything about the quality of the given structure, but are here used as a way to avoid issues when automating simulations. 
                   



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
}

