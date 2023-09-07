import MDAnalysis as mda
from numpy.linalg import inv
from numpy import array
from numpy import matrix
from numpy import float32
import sys
import argparse

##
## you need the following libraries:
##  conda install mdanalysis -c conda-forge
##
## Written by Max Bonomi (mbonomi@pasteur.fr)
##

# initialize parser
parser = argparse.ArgumentParser(prog='python align-XTC.py', description='Align XTC trajectory based on a pre-determined transformation')
parser.add_argument('in_pdb',    type=str, help='input PDB file')
parser.add_argument('in_xtc',    type=str, help='input XTC trajectory to align')
parser.add_argument('out_xtc',   type=str, help='output XTC trajectory file')
parser.add_argument('transform', type=str, help='file with transformation to apply')
args = parser.parse_args()

#### INPUT
# input PDB
IN_PDB=vars(args)["in_pdb"]
# input XTC 
IN_XTC=vars(args)["in_xtc"]
# output XTC
OUT_XTC=vars(args)["out_xtc"]
# file with transformation
TRASF=vars(args)["transform"]

# read transformation from file (in a dictionary - Angstrom)
transf={}
f = open(TRASF, "r")
transf = eval(f.read())

# calculate inverse of transformation matrix
Rinv = inv(transf["R"]) 

# create universe
mobile = mda.Universe(IN_PDB,IN_XTC)
# select all atoms
atoms = mobile.select_atoms('all')

# new trajectory
W = mda.Writer(OUT_XTC, atoms.n_atoms)

# cycle on frames
for i in range(0, len(mobile.trajectory)):
   # go to frame i
   mobile.trajectory[i]
   # transform atoms
   atoms.translate(-transf["X1"])
   atoms.rotate(Rinv)
   atoms.translate(transf["X0"])
   # write frame
   W.write(atoms)
