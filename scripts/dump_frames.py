import os,glob, json
import mdtraj as md
import argparse


def dump_traj(traj_fn, pdb_fn, freq = None, o_dir = './traj_pdbs' ):
    '''
    Takes trajectory file and cleaned pdb file and splits
    trajectory into n frames that are each save as a pdb file
    '''
    assert os.path.exists(o_dir), 'Output directory does not exit'
    
    name_pdb = pdb_fn.split('/')[-1].split('.pdb')[0]
    
    traj = md.load(traj_fn, top=pdb_fn)
    
    # of no timestep freq, output all frames in traj
    if freq == None: 
        for frame in traj:
            t = int(frame.time[0])
            o_pdb = o_dir + f'/{name_pdb}_{t}ps.pdb'
            frame.save_pdb(o_pdb) 
            
    # output frames every freq step
    else:
        for frame in traj:
            t = int(frame.time[0])
            if int(t) % freq == 0:
                o_pdb = o_dir + f'/{name_pdb}_{t}ps.pdb'
                frame.save_pdb(o_pdb) 
    
        # save last frame in traj
        last = int(traj.time[-1])
        o_pdb = o_dir + f'/{name_pdb}_{last}ps.pdb'
        frame.save_pdb(o_pdb) 


def parser():
    '''
    makes parser for script
    '''
    
    parser = argparse.ArgumentParser(description='Converts trajectory to single framed pdbs')
    parser.add_argument('-traj',
                        dest='traj',
                        type=str,  
                        help='name of traj file')
    parser.add_argument('-pdb',
                        dest='pdb',
                        type=str,
                        help='path to pdb topology file, (only protein no solvent)'
                       )
    parser.add_argument('-freq',
                        dest='freq',
                        default = None,
                        type=int,
                        help='how often to output frames from simulation in ps'
                       )
    parser.add_argument('-o_dir',
                        dest='o_dir', 
                        type=str,
                        default='./traj_pdbs',
                        help='output directory (default: traj_pdbs)')

    
    args = parser.parse_args()

    return args


    
if __name__=='__main__':
    
    args = parser()
    
    print(f'Splitting trajectory {args.traj}. Dumping frames every {args.freq} ')
    
    # create output directory
    if not os.path.isdir(args.o_dir):
        os.mkdir(args.o_dir)
    print(f'output pdb frames to be found in {args.o_dir}')
    
    dump_traj(traj_fn = args.traj, 
               pdb_fn = args.pdb,
               freq = args.freq,
               o_dir = args.o_dir)
    
    
    
