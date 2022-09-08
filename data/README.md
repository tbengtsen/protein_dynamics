# DATA INFORMATION
The PDB data used for these simulations are based on the CATH 4.2 [1] Protein Structure Classification database which filters and classifies protein domains based on evolutionary relationship. This is used to obtain a balanced dataset. The dataset split used here are copied from John Ingraham's paper [2] for later comparison of results. 

In addition, we here filter out proteins from the CATH 4.2 based on whether their structure are appropiate for these types of automated MD simulations. 

# Furtheer data filter criterias for simulations
Not all PDB files present in the CATH 4.2 database are suitable for these kinds of automated molecular dynamics simulations. We therefor apply further filtering criteria as described below:

    - if it is a membrane protein => excluded
    - if it is an unsuitable experimental method  for solving the structure => excluded (unsuitable only in the sense of inappropiate for these automated simulations). Here the following methods are included and all others hence excluded:
        - ELECTRON MICROSCOPY => included
        - SOLID-STATE NMR => included
        - SOLUTION NMR => included
        - X-RAY DIFFRACTION => included
        - ELECTRON CRYSTALLOGRAPHY => included
    - if the X-RAY resolution is >=2.5  => excluded
    - to high RMSD (>3.5) between models in NMR models => excluded 
    - contains other amino acids  than the standard 20 amino acids and Selenocystein => excluded
    - gaps/unresolved parts inside the chain of a protein (in terminals these are allowed) => excluded
    - unsuitable ligand(s) (unsuitable only in the sense of usefull for automation of simulations) such as NATP => excluded
    - too big a protein i.e >500 residues. This threshold keeps the computational ressources lower for simulations and keeps memory usage down in ML algorithms => excluded 
    - any chain too small in the  protein (<50 residues). Meaning it is unlikely to process a (somewhat) 3D structure by itself => excluded

**OBS!** It should be stressed that neither of the filtering conditions descriped above necessarily say anything about the quality of the given structure, but are here used as a way to avoid issues when automating simulations. 
                   
For more details, see the filter script in `CATH4.2/filter_cath_dataset.ipynb`

# PDB data processing
## Download format
The downloaded PDB files are downloaded in their bioassembly form from `ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates` if existing.

## Rebuild missing atoms or residueese
As many PDB files have missing atoms or even residues these are rebuildt into the PDB file. NB as mentioned above only terminal missing residues are rebuilt here. However, a script exists to rebuild gaps in the structure `rebuild_by_modeller.py` in the `scripts/clean_pdb_functions/` directory. This is not automated and thus only used for downstream calculations of stability data.

The programs Reduce [3] together with pdbfixer from the great OpenMM [4] are used to estimate charges, and for rebuilding other atoms and residues. For more information, see the scripts in scripts/clean_pdb_functions/` directory







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

3)
@article{Reduce,
author = {Word, et. al.},
title = {Asparagine and Glutamine: Using Hydrogen Atom Contacts in the Choice of Side-chain Amide Orientation},
journal = { J. Mol. Biol.}
year = {1999}, 
volumen = {285 },
pages = { 1733-1747},
website = {https://github.com/rlabduke/reduce}}

4)
@article{OpenMM,
author = {P. Eastman, J. Swails, J. D. Chodera, R. T. McGibbon, Y. Zhao, K. A. Beauchamp, L.-P. Wang, A. C. Simmonett, M. P. Harrigan, C. D. Stern, R. P. Wiewiora, B. R. Brooks, and V. S. Pande. },
titleee = { OpenMM 7: Rapid development of high performance algorithms for molecular dynamics},
journal = {PLOS Comp. Biol.},
year = {2017},
volumen = { 13(7): e1005659. },
website ={ https://github.com/openmm/pdbfixer }}
