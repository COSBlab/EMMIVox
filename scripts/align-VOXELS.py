import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import numpy as np
import sys
import argparse

##
## you need the following libraries:
##  conda install mdanalysis -c conda-forge 
##
## Written by Max Bonomi (mbonomi@pasteur.fr) and Samuel Hoff (shoff@pasteur.fr)
##

# initialize parser
parser = argparse.ArgumentParser(prog='python align-VOXELS.py', description='Align mobile to reference PDB and transform the map associated to the mobile PDB')
parser.add_argument('ref_pdb',      type=str, help='MD PDB file (reference)')
parser.add_argument('mobile_pdb',   type=str, help='RCSB PDB file to align to ref_pdb (mobile)')
parser.add_argument('ref_map',      type=str, help='PLUMED output map aligned to ref_pdb')
parser.add_argument('mobile_map',   type=str, help='PLUMED input map aligned to mobile_pdb')
parser.add_argument('--ref_sel',    type=str, nargs=1, help="MDAnalysis selection for alignment: MD PDB atoms")
parser.add_argument('--mobile_sel', type=str, nargs=1, help="MDAnalysis selection for alignment: RCSB PDB atoms")
args = parser.parse_args()

#### INPUT
# reference PDB
REF_PDB=vars(args)["ref_pdb"]
# mobile PDB to align
MOB_PDB=vars(args)["mobile_pdb"]
# output map
OUT_MAP=vars(args)["ref_map"]
# input map to transform
IN_MAP=vars(args)["mobile_map"]
# selection
selection='backbone'; selection2='backbone'
if(vars(args)["ref_sel"]!=None):
  selection=vars(args)["ref_sel"][0]
if(vars(args)["mobile_sel"]!=None):
  selection2=vars(args)["mobile_sel"][0]

# create reference and mobile universes
ref    = mda.Universe(REF_PDB)
mobile = mda.Universe(MOB_PDB)
# define two sets of atoms based on selection
ref_atoms    = ref.select_atoms(selection)
mobile_atoms = mobile.select_atoms(selection2)
# print to output to help with debugging
print("Reference atom count: ", len(ref_atoms))
print("Mobile atom count: ",    len(mobile_atoms))
# check lengths
if(len(ref_atoms) != len(mobile_atoms)):
  print("ERROR: Reference and mobile selections have a different number of atoms")
  exit()
# find transformation that minimizes RMSD
ref0    =    ref_atoms.positions -    ref_atoms.center_of_geometry()
mobile0 = mobile_atoms.positions - mobile_atoms.center_of_geometry()
R, rmsd = align.rotation_matrix(mobile0, ref0)
print("RMSD [Ang]: ", rmsd)

# store transformation to file (in a dictionary - Angstrom)
transf={}
transf["X0"] = mobile_atoms.center_of_geometry()
transf["R"]  = R
transf["X1"] = ref_atoms.center_of_geometry()
f = open("transformation.dat", "w")
f.write(str(transf)); f.close()

# open input MAP
pos=[]; density=[]; error=[]
for lines in open(IN_MAP, "r").readlines():
    riga=lines.strip().split()
    # check header
    if(riga[0]=="#!"): 
       header=lines.strip()
       continue
    # convert coordinates to Ang
    x = 10.0 * float(riga[1])
    y = 10.0 * float(riga[2])
    z = 10.0 * float(riga[3])
    # add coordinates to list
    pos.append([x, y, z])
    # add density to list
    density.append(float(riga[4]))
    # add error to list
    error.append(float(riga[5]))

# transfor pos to numpy array
pos = np.array(pos)
# translate COM of mobile
pos -= transf["X0"] 
# rotate
pos = np.transpose(np.matmul(transf["R"],np.transpose(pos)))
# translate COM of ref 
pos += transf["X1"]

# print new MAP
omap=open(OUT_MAP, "w")
# write header
omap.write("%s\n" % header)
# cycle on voxels
for i in range(0, len(density)):
    # print out
    omap.write("          %6d  " % i)
    for j in range(0,3):
        # back to nanometers first
        x = 0.1 * pos[i,j]
        omap.write("%13.7lf " % x)
    omap.write("  %10.7e  %10.7e\n"  % (density[i],error[i]))
# close file
omap.close()
