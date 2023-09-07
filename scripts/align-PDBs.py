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
parser = argparse.ArgumentParser(prog='python align-PDBs.py', description='Align PDBs based on a pre-determined transformation')
parser.add_argument('in_pdb',    type=str, help='input PDB file to align')
parser.add_argument('out_pdb',   type=str, help='output PDB file')
parser.add_argument('transform', type=str, help='file with transformation to apply')
args = parser.parse_args()

#### INPUT
# input PDB
IN_PDB=vars(args)["in_pdb"]
# output PDB 
OUT_PDB=vars(args)["out_pdb"]
# file with transformation
TRASF=vars(args)["transform"]

# read transformation from file (in a dictionary - Angstrom)
transf={}
f = open(TRASF, "r")
transf = eval(f.read())

# calculate inverse of transformation matrix
Rinv = inv(transf["R"]) 

# create universe
mobile = mda.Universe(IN_PDB)

# transform atoms and write PDB 
mobile.atoms.translate(-transf["X1"])
mobile.atoms.rotate(Rinv)
mobile.atoms.translate(transf["X0"])
mobile.atoms.write(OUT_PDB)
