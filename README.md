# Automated  protein dynamics protocol 

This repo contains the means to automate short molecular dynamics simulations on protein structures from a raw PDB file to an equilibrated and then simulated trajectory output file. It is based on the great OpenMM simulation package [[1](#openmm)]. 

The repo was created to obtain a dynamics augmented protein structure dataset to test the impact of using dynamics augmentation as input to train  machine learning models and to test the impact of dynamics for improvement in downstream protein stability predictions.  The 1282 finished 20 ns simulations and their outputs are currently not freely available but this is being worked on. For access, contact Wouter Boomsma at University of Copenhagen, and for more information about the dataset see the below in:  [Augmentated Dynamics Dataset](dataset_info)

The whole protocol is contained in the script `protocol_molecular_dynamics_simulations.py` or as a notebook `scripts/notebook_protocol/openMM_protocol.ipynb` with helper functions in the scripts directory. See below for necessary [installations](installations), [environment](environment) and [how to run](how_to_run). 


# Automated protocol 
The protocol automates the simulations from a start pdb structure to end output simulation trajectory. It can be run as a seperate command line python script or in a notebook and performs the following steps:

    1) Preps/cleans a PDB file and rebuild missing atoms and residues to set up for simulation (obs not part of protocol script but seperate script used)
        
    2) Builds the simulation system with 3 choices of simulation box, addition of water and ions for neutralising system
    
    3) Equilibrates the system in 3 steps:
        a) Minimization w. restraints on heavy atoms in backbone and sidechains 
        b) 1st equilibration with restraints  on heavy atoms in backbone and sidechains and a rolling window to determine whether equilibration is reached
        c) 2nd equilibration without restraints  on any atoms a rolling window to determine whether equilibration is reached

    4) Performs the production run (actual simulation) for 20 ns (default) with an inbuildt check to see if simulation blew up 
    
    5) Outputs simulation trajectory with a choice of output file (default XTC without water)
    
    6) View single simulation in notebook using the awesome NGLview package [4](see  `scripts/notebook_protocol/view_single_simulation.ipynb` )


For more information on the detailed specifics of the simulation msuch as force field model, restraint size etc, please see below: [Methods - details on the simulation](simulations_details)  




# Installations
<a name="installations"></a>
Install the program Reduce [3] seperately. This is used to clean PDB files and is used by the parser to set charges and add missing hydrogens to the proteins.
```
cd scripts/
git clone https://github.com/rlabduke/reduce.git 
cd reduce/
make; make install # This might give an error but provide the reduce executable in this directory
```

# Environment
<a name="environment"></a>
See the file: `environment.yml` for which python packages and versions is used. 




# How to run
<a name="how_to_run"></a>
1) First inspect the PDB file to see if there is larger gaps in the structure. If so, it is recommended to use the Modeller [5]clean the PDB file using the script : 




# Methods - details on the simulation 
<a name="simulations_details"></a>
The details of each step in the automated protocol are as follows. Some steps have multiple options, some does not for various reasons. 


## 1)  Prep PDB input files
Preps downloaded PDB files to be able to parse to simulation protocol. It cleans most ligands from the structure and rebuilds missing atoms and, if specified, rebuilds missing amino acids. 
This part is seperated from the actual simulation protocol. The scripts needed are found in the: `scripts/clean_pdb_functions/`. 


### Download the correct conformation
Firstly, it is strongly recommended to download a bioassembly PDB file, i.e. a PDB that only contains the actual naturally occuring protein structure and not e.g. X-ray asymmetrical unit. See more at [RSCB's guide to biological assemblies](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/biological-assemblies)

### Clean PDB file
Cleans the downloaded PDB file(s) to be able to parse to simulation protocol. This includes a multitude of hacks that each decrease the number of time the parsing crashes when running on a large number of pdb files. Not all steps are necessary in most cases, and often only the Non heterogenselector, the pdbfixer and the last reduce step is necessary.

The is three scripts available dependent in the needs: 
    1) do not rebuild any missing residues but clean and rebuild missing atoms:
        - use `clean_without_rebuild.py`
    2) rebuild missing residues with Modeller (needs some manual work but might give lower chances of crashes in the simulation)
        - a) use the `rebuild_by_modeller.py` and remember to make alignment file manually
        - b) use the output of a) to clean and rebuild missing atoms in `clean_without_rebuild.py`
    3) rebuild missing residues with pdbfixer and clean and rebuilds missing atoms (observed 
        - use the `clean_and_rebuild_by_pdbfixer.py` 
        
Here it is recommended to use the more cumbersome Modeller based approach as the pdbfixer based have (personal observations only) been observed to give som strange conformations of e.g. Prolines, Alanines and Threonines in the rebuilt sections. However, this might work with proper minimization, but also just crashing the simulation during minimization. 


## 2)  Build simulation system
This steps builds the cleaned protein structure into a simulation box with surrounding water and ions for neutralising system, and a padding (distance from protein to box edge). 
The choice of geometry of the box is important for the computational ressources needed as the smaller the volumen, the fewer atoms and, hence, the fewer calculations of forces in each timestep as well as how well the system scales the calculations. In addition, it is important if the protein is not globular and the simulation is long enough for the protein to rotate - here, a non-symmetric box should not be chosen. 
s
### Default system setups

    - forcefield: Amber [ff14SB](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00255) 
    - water type: TIP3P 
    - box geometry : rhombic dodecahedron
    - padding : 1 nm on each side of solvated protein i.e. 2 nm in total
    
For more information on installation of forcefields and use of forcefields see [OpenMM](https://github.com/openmm/openmmforcefields)
 
For more information on the math behind how trigonal boxes can be simulated as cubic to be computional efficient see  [paper below](box_geometry). 


## 3) Equlibrate system


### Default equilibration setup

# Augmented dynamics dataset
<a name="dataset_info"></a>

The input dataset of proteins where obtained from CATH 4.2  [2] Protein Structure Classification database which filters and classifies protein domains based on evolutionary relationship. This was used to obtain a balanced dataset.

The automated simulation data contains:
    - 2195 of the PDBs in CATH4.2 that was deemed appropiate for simulations based on filtering criteria descriped in the README.md file in data/
    - 1282 finished simulations of proteins from CATH4.2 (out of the 2195 attempted simulations)
    - of which 505 PDBs are obtained from NMR structures
    - of which 777 PDBs are obtained from X-ray structures
    - of the 911 PDBs that was not simulated: 
        - 522 was not included due to too high RMSD at the end of the simulation
        - 161 PDBs was not equilibrated during the equilibration steps 
        - 159 PDBs blew up with NaN errors during equilibration or simulation, typically because they were not appropiate for simulations but not cought by our filter criteria
        - 69 had problems with I/O that persisted to the end of this postdoc
    





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



# Citations

<a name="openmm"></a>

```
1) 
@article{OpenMM,
author = {P. Eastman, J. Swails, J. D. Chodera, R. T. McGibbon, Y. Zhao, K. A. Beauchamp, L.-P. Wang, A. C. Simmonett, M. P. Harrigan, C. D. Stern, R. P. Wiewiora, B. R. Brooks, and V. S. Pande. },
title = { OpenMM 7: Rapid development of high performance algorithms for molecular dynamics},
journal = {PLOS Comp. Biol.},
year = {2017},
volumen = { 13(7): e1005659. },
website ={ https://github.com/openmm/pdbfixer }}
```

```2) 
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
3)
@article{Reduce,
author = {Word, et. al.},
title = {Asparagine and Glutamine: Using Hydrogen Atom Contacts in the Choice of Side-chain Amide Orientation},
journal = { J. Mol. Biol.}
year = {1999}, 
volumen = {285 },
pages = { 1733-1747},
website = {https://github.com/rlabduke/reduce}}
```


```
4)
@article{NGLview,
    author = {Nguyen, Hai and Case, David A and Rose, Alexander S},
    title = "{NGLview–interactive molecular graphics for Jupyter notebooks}",
    journal = {Bioinformatics},
    year = {2017},
    volume = {34},
    number = {7},
    pages = {1241-1242},
    doi = {10.1093/bioinformatics/btx789},
    website = {https://github.com/nglviewer/nglview}
}


```

```
5) 
@article{Modeller,
author = {Webb, Benjamin and Sali, Andrej},
title = {Comparative Protein Structure Modeling Using MODELLER},
journal = {Current Protocols in Bioinformatics},
year = {2014}
volume = {47},
number = {1},
pages = {5.6.1-5.6.32},
doi = {https://doi.org/10.1002/0471250953.bi0506s47}}

```

<a name="box_geometry"></a>
```
@article{simulation_box,
author = {Bekker, H.},
title = {Unification of box shapes in molecular simulations},
journal = {Journal of Computational Chemistry},
year = {1997}
volume = {18},
number = {15},
pages = {1930-1942},
doi = {https://doi.org/10.1002/(SICI)1096-987X(19971130)18:15<1930::AID-JCC8>3.0.CO;2-P},
}

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


