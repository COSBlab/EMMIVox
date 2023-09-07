import MDAnalysis as mda
import numpy as np
import sys
import argparse

##
## you need the following libraries:
##  conda install mdanalysis -c conda-forge
##
## Written by Max Bonomi (mbonomi@pasteur.fr)
##

# initialize parser
parser = argparse.ArgumentParser(prog='python add-BFACT.py', description='Write PDB with Bfactor extracted from PLUMED status file')
parser.add_argument('in_pdb',  type=str, help='input PDB file')
parser.add_argument('bfactor', type=str, help='PLUMED status file with Bfactor')
parser.add_argument('out_pdb', type=str, help='output PDB file')
parser.add_argument('--exclude',   type=str, nargs=1, help="MDAnalysis selection of atoms to exclude in output PDB")
args = parser.parse_args()

#### INPUT
# input PDB
IN_PDB=vars(args)["in_pdb"]
# bfactors
BFACT=vars(args)["bfactor"]
# output PDB
OUT_PDB=vars(args)["out_pdb"]
# selection
SEL="all"
if(vars(args)["exclude"]!=None):
  SEL += " and not ( "+vars(args)["exclude"][0] +" )"

# create universe
ref=mda.Universe(IN_PDB)
# select atoms
sel=ref.select_atoms(SEL)

# read PLUMED bfactor file
res_c={}; bfact={}; nochain=False; found=False
for lines in open(BFACT, "r").readlines():
    riga=lines.strip().split()
    # line with legend
    if(riga[0]=="#!"):
      # store keys
      for i in range(0,len(riga)):
          # check if "bf-" in riga[i]
          if("bf-" not in riga[i]):
            continue
          else:
            if(found==False):
              # store column where bfactors start
              j=i-2
              # set flag
              found=True
          # get residue,chain pair
          key=riga[i].split("-")[1]
          # residue id
          res=int(key.split(":")[0])
          # chain
          c=key.split(":")[1]
          # check if nochain
          if(c=="*"): nochain=True
          # map column pair res,c
          res_c[i-2] = (res,c)
    # otherwise, store bfactor dictionary
    # last line of status file will be retained
    else:
      for i in range(j,len(riga)):
          # convert unit of measures
          bfact[res_c[i]] = 100.0 * float(riga[i])

# now modify bfactor field
for at in sel:
    # get chain
    c = at.segid
    if(nochain): c="*"
    # set temp factor
    at.tempfactor = bfact[(at.resid,c)]

# and save pdb file
sel.write(OUT_PDB)
