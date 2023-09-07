import MDAnalysis as mda
import math
import numpy as np
import argparse

##
## you need the following libraries:
##  conda install mdanalysis -c conda-forge
##
## Written by Max Bonomi (mbonomi@pasteur.fr) and Samuel Hoff (shoff@pasteur.fr)
##

# initialize parser
parser = argparse.ArgumentParser(prog='python make_ndx.py', description='Create GROMACS index.ndx file with atoms selection')
parser.add_argument('structure',   type=str, help='PDB/GRO file')
parser.add_argument('selection',   type=str, help='MDAnalysis selection')
parser.add_argument('group_name',  type=str, help='Group name')
parser.add_argument('--ndx',       type=str, nargs=1, help='GROMACS ndx file to append output')
parser.add_argument('--water',     action='store_true', help='Swap water positions to have selection at the top of file')
args = parser.parse_args()

#### INPUT
# input PDB/GRO file
PDB = vars(args)["structure"]
# selection
selection = vars(args)["selection"]
# group name
gname = vars(args)["group_name"] 
# output file name
oname = "index.ndx"
if(vars(args)["ndx"]!=None):
  oname=vars(args)["ndx"][0]
# swap waters
water = vars(args)["water"]

# load PDB
u = mda.Universe(PDB)
# atom selection - exclude hydrogens and COO-
atoms = u.select_atoms("( "+selection+' ) and not type H and not (resname GLU ASP and name OD1 OD2 OE1 OE2)', periodic = False)

# check number of waters
resn = [r.resname for r in atoms.residues]
nwat = resn.count('TIP3')+resn.count('TIP4')+resn.count('TIP5')+resn.count('SOL')+resn.count('HOH') 
if(nwat>0 and nwat!=len(resn) and water):
  print('ERROR: you should not specify the flag --water if you selected both water and non-water atoms')
  exit()
if(nwat==len(resn) and not water):
  print('ERROR: you selected only water atoms but did not add the flag --water')
  exit()

# swap waters - hydrogen included
if(water):
  # get all atoms
  atoms = u.atoms
  # all waters
  wat = atoms.select_atoms('resname '+resn[0])
  # cryo-EM waters
  wat0 = atoms.select_atoms(selection, periodic = False)
  # other waters
  wat1 = wat - wat0
  # do the swap
  j=0
  for i in range(len(wat0)):
      # check if you need to swap
      if(wat0[i].resid < wat1[j].resid): continue
      pos = wat1[j].position
      wat1[j].position = wat0[i].position
      wat0[i].position = pos
      # increase j
      j+=1
  # write final gro file
  atoms.write(PDB)
  # re-load PDB
  u = mda.Universe(PDB)
  # redo atom selection - exclude hydrogens and COO-
  atoms = u.select_atoms("( "+selection+' ) and not type H and not (resname GLU ASP and name OD1 OD2 OE1 OE2)', periodic = False)

# append to index file
out=open(oname, "a")
out.write("[ "+gname+" ]\n")
for i,at in enumerate(atoms):
    out.write("%7d " % (at.index+1))
    # go to new line
    if((i+1)%15==0 or i==len(atoms)-1): out.write("\n")

# add a second group with original selection (with hydrogen atoms and COO-)
atoms = u.select_atoms(selection, periodic = False)
# write to index
out.write("[ "+gname+"-H ]\n")
for i,at in enumerate(atoms):
  out.write("%7d " % (at.index+1))
  # go to new line
  if((i+1)%15==0 or i==len(atoms)-1): out.write("\n")
