# Different methods to clean a pdb 
The biggest difference is whether to rebuild missing residues in a pdb 
AND how to do so. 
Here it is recommended to use the more cumbersome Modeller based approach
as the pdbfixer based have here been observed to give som strange 
conformations of e.g. Prolines, Alanines and Threonines in the rebuilt sections
However, this might work with proper minimization, but also just risk 
crashing the simulation during minimization 

## functions to clean PDBs

    1) do not rebuild any missing residues
        - use `clean_without_rebuild.py`
    2) rebuild missing residues with Modeller
        - a) use the `rebuild_by_modeller.py` and remember to make alignment file manually
        - b) use the output of a) in `clean_without_rebuild.py`
    3) rebuild missing residues with pdbfixer 
        - use the `clean_and_rebuild_by_pdbfixer.py` 
        