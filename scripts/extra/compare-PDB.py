import MDAnalysis as mda
import numpy as np
import sys
import argparse
import math

##
## you need the following libraries:
##  conda install mdanalysis -c conda-forge
##
## Written by Max Bonomi (mbonomi@pasteur.fr)
##

def read_cross_corr(ccfile):
    # initialize dictionary
    cc={}
    # read file
    for lines in open(ccfile, "r").readlines():
        riga=lines.strip().split()
        # define key: (residue, chain) 
        key = (int(riga[2]),riga[0])
        # add to dictionary
        cc[key] = float(riga[3])
    return cc

# initialize parser
parser = argparse.ArgumentParser(prog='python compare-PDB.py', description='Compare two PDBs - cc per residue and RMSD')
parser.add_argument('in_pdb',  type=str, help='input PDB file')
parser.add_argument('cc',      type=str, help='cc per residue of input PDB')
parser.add_argument('out_pdb', type=str, help='output PDB file')
parser.add_argument('--ccref', type=str, help='cc per residue of reference model')
parser.add_argument('--sel',   type=str, nargs=1, help="MDAnalysis selection for output PDB")
args = parser.parse_args()

#### INPUT
# input PDB
IN_PDB=vars(args)["in_pdb"]
# cross-correlation per residue of input PDB
CC=vars(args)["cc"]
# output PDB
OUT_PDB=vars(args)["out_pdb"]
# cross-correlation per residue of reference model
CCREF=vars(args)["ccref"]
# compare with reference cross correlation
do_compare=False
if(CCREF!=None): do_compare=True
# selection
SEL="not type H"
if(vars(args)["sel"]!=None):
  SEL=vars(args)["sel"][0]

# create universe
pdb=mda.Universe(IN_PDB)
# select atoms
sel=pdb.select_atoms(SEL)
# get dictionary cross-correlation per residue of input PDB
cc_dict = read_cross_corr(CC)
# get dictionary cross-correlation per residue of reference model
# and compare
if(do_compare): 
  cc_dict_ref = read_cross_corr(CCREF)
  # calculate difference for keys in cc_dict
  for key in cc_dict:
      # check if in cc_dict_ref
      if key in cc_dict_ref:
         cc_dict[key] -= cc_dict_ref[key]
      # otherwise set to zero
      else:   
         cc_dict[key] = 0.0

# now modify bfactor field
for at in sel:
    # set temp factor
    at.tempfactor = cc_dict[(at.resid,at.segid)]

# and save pdb file
sel.write(OUT_PDB)
