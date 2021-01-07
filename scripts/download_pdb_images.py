#!/usr/bin/env python 

'''
script to download pictures of pdbs from rcsb to check why failed
'''

import requests
import sys
import os
import argparse

def parser():

    parser = argparse.ArgumentParser(description='download pdb images from rcsb.')
    parser.add_argument('-list',
                        dest='list_pdbs',
                        type=str,  
                        help='iname of file with list of pdbs id, only 4 letter id')
    parser.add_argument('-o_dir',
                        dest='o_dir', 
                        type=str,
                        default='images_pdbs',
                        help='output directory (default: images_pdbs)')

    args = parser.parse_args()

    return args

def read_pdbs(fn):

    with open(fn, 'r') as f:
        pdbs = f.readlines()
    pdbs = [pdb.strip() for pdb in pdbs]
    # check only pdb ids:
    for pdb in pdbs:
        assert len(list(pdb)) == 4, 'obs input list file must only contain pdb ids of 4 letters'
    return pdbs

def download_pictures(pdbs, o_dir):
    'downloads all pdb pictures from rcsb to output directory'
    if not os.path.isdir(o_dir):
        os.mkdir(o_dir)
    for pdb in pdbs:
        mid_letters = list(pdb)[1:3]
        mid_letters = ''.join (mid_letters)
        image_url = f'https://cdn.rcsb.org/images/structures/{mid_letters}/{pdb}/{pdb}_assembly-1.jpeg'
        img_data = requests.get(image_url).content
        with open(f'{o_dir}/image_{pdb}.jpg', 'wb') as handler:
                handler.write(img_data)
    
     
if __name__=="__main__":
    
    args = parser()
    pdbs = read_pdbs(args.list_pdbs)
    download_pictures(pdbs, args.o_dir)



