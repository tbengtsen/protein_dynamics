'''
Example script for running Sali's Modeller to rebuilt missing residues 
in a pdb structure file. 
This demands an alignment file where its detailed the missing residues. 
See how to build this at: 
https://salilab.org/modeller/8v2/manual/node176.html

OBS! This rebuild method demands a license from modeller, see:
https://salilab.org/modeller/registration.html
Its free for academics. 
In addition it demands the following:
    - download or use conda to install the Modeller package

'''

import argparse, os
from modeller import log
from modeller import environ
from modeller import alignment
from modeller.automodel import *


if __name__ == "__main__":
    # Argument Parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_file_in", type=str)
    parser.add_argument("--alignment_file",
                        type=str,
                        required=True,
                        help='modellers alignment file. Must be created manually.') 
    parser.add_argument('--clean_up', 
                        action="store_true", 
                        help='whether to remove modellers non pdb output files')
    
    # Parse arguments
    args_dict = vars(parser.parse_args())
    pdb_input_filename = args_dict["pdb_file_in"]
    alignment_file = args_dict["alignment_file"]
    clean_up = args_dict["clean_up"]
    
    # set output name 
    o_name  =  pdb_input_filename.split('.pdb')[0]
    o_name = o_name + '_rebuilt.pdb'
 
    # Setup modeller
    log.verbose()
    env = environ()
    env.io.atom_files_directory = ['.', '.']
    a = automodel(env,
                  alnfile  = alignment_file,           # alignment filename
                  knowns   = 'STRUCT',              # name of the template structure set in align file
                  sequence = 'SEQ')               # name of entire sequence defined in align file
    
    a.starting_model= 1                 # index of the first model
    a.ending_model  = 1                 # index of the last model
                                        # (determines how many models to calculate)
    a.make()
    a.write(file=o_name, model_format='PDB')

    if clean_up:
        os.system('rm *.ini *.rsr *.sch  *.V99990001 *.B99990001.pdb *.D00000001')
