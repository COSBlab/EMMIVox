import MDAnalysis as mda
import math
import numpy as np
import argparse
from scipy.spatial import distance

##
## you need the following libraries:
##  conda install mdanalysis -c conda-forge
##
## Written by Max Bonomi (mbonomi@pasteur.fr) and Samuel Hoff (shoff@pasteur.fr)
##

# get atoms indexes from index file
def get_atom_indexes(oname,gname):
    idx=[]; toread=False
    for line in open(oname, "r").readlines():
        riga=line.strip().split()
        # empty line
        if(len(riga)==0): continue
        # detect start of the group
        if(riga[0]=="[" and riga[1]==gname): toread=True
        if(riga[0]=="[" and riga[1]!=gname): toread=False 
        # read line if right group
        if(riga[0]!="[" and toread):
           for i in range(0, len(riga)):
               idx.append(int(riga[i])-1)
    return idx 

# write atoms indexes to file
def write_to_ndx(ifile, gname, idx):
    out = open(ifile, 'a')
    # write to index
    out.write("[ "+gname+" ]\n")
    for i,j in enumerate(idx):
        out.write("%7d " % j)
        # go to new line
        if((i+1)%15==0 or i==len(idx)-1): out.write("\n")
    out.close()

# initialize parser
parser = argparse.ArgumentParser(prog='python add-POSRES.py', description='Add position restraint to water molecules')
parser.add_argument('structure',   type=str,   help='PDB/GRO file')
parser.add_argument('selection',   type=str,   help='MDAnalysis selection for ordered water molecules')
parser.add_argument('upper_limit', type=float, help='Distance upper limit for water oxygens')
parser.add_argument('--ndx',       type=str, default="index.ndx", help='GROMACS ndx file')
parser.add_argument('--equil',     action='store_true', help='Positional restraints for equilibration')
parser.add_argument('--kappa',     type=float, default=10000, help='GROMACS ndx file')
args = parser.parse_args()

#### INPUT
# input PDB/GRO file
PDB = vars(args)["structure"]
# selection of ordered water molecules
selection = vars(args)["selection"]
# upper limit distance (convert to nm)
dmax = vars(args)["upper_limit"] * 0.1
# output file name
oname=vars(args)["ndx"]
# flag for equilibration
doEquil = vars(args)["equil"]
# set stride
stride=4
if(doEquil): stride=1
# restraint intensities 
KAPPA_=vars(args)["kappa"]

# load PDB
u = mda.Universe(PDB)
# get ordered water molecules (exclude hydrogens)
# from MDAnalysis selection
atoms_ow = u.select_atoms(selection+" and not type H")
# get indexes
idx_ow = [at.index for at in atoms_ow]
# get water atoms indexes from file: ordered molecules plus buffer
idx_w = get_atom_indexes(oname,"System-WAT")
# and the corresponding atoms
atoms_w = u.atoms[np.array(idx_w)]
# get protein atoms indexes from file
idx_m = get_atom_indexes(oname,"System-PRO")
# and the corresponding atoms (select backbone - more rigid during equilibration)
atoms_m = u.atoms[np.array(idx_m)].select_atoms('backbone or nucleicbackbone')

# now for each atom in atoms_w, find closest atoms_ow/atoms_m
w_close={}; ow_dict={}
for at in atoms_w:
    # distances between atom and all ordered waters 
    dist=distance.cdist([at.position], atoms_ow.positions).squeeze()
    # find closest ordered water
    w = atoms_ow[np.argmin(dist)]
    # distances between ordered water and all protein atoms 
    dist=distance.cdist([w.position], atoms_m.positions).squeeze()
    # if this is an ordered water and equilibration
    if(at.index in idx_ow and doEquil):
      # order by distance 
      ic = np.argsort(dist)
      # add 3 points at different distances to dictionary
      for i in [0, int(len(ic)/2), len(ic)-1]:
          # protein atom
          m = atoms_m[ic[i]]
          # distance from protein atom
          d = np.linalg.norm(at.position-m.position)
          # convert to nm and store in dictionary
          ow_dict[(at.index+1,m.index+1)] = d * 0.1 
    # find the closest protein atom for buffer restraint
    m = atoms_m[np.argmin(dist)]
    # store in dictionary
    w_close[at.index+1] = m.index+1

# create plumed posres file
out=open("plumed_posres.dat", "w")
# Equilibration 
if(doEquil):
   # create groups for WRAPAROUND action
   g0=[]; g1=[]; g2=[]
   # loop over all ordered molecules plus buffer
   for key in w_close:
       # ATOMS
       g0.append(key)
       # with hydrogens
       g1.append(key); g1.append(key+1); g1.append(key+2)
       # AROUND
       g2.append(w_close[key])
   # write groups to index-wrap file
   write_to_ndx("index-wrap.ndx", "Wrap-ATOMS",   g0)
   write_to_ndx("index-wrap.ndx", "Wrap-ATOMS-H", g1)
   write_to_ndx("index-wrap.ndx", "Wrap-AROUND",  g2)
   # add header
   out.write("# include topology info\n")
   out.write("MOLINFO STRUCTURE=../3-Map-Scaling/step3_input_xtc.pdb WHOLE\n")
   out.write("# define protein atoms\n");
   out.write("system-pro: GROUP NDX_FILE=../0-Building/index.ndx NDX_GROUP=System-PRO\n")
   out.write("# make protein whole\n");
   out.write("WHOLEMOLECULES ADDREFERENCE EMST ENTITY0=system-pro\n")
   out.write("# Wrap water around protein\n")
   out.write("wrap-at: GROUP NDX_FILE=index-wrap.ndx NDX_GROUP=Wrap-ATOMS\n")
   out.write("wrap-ar: GROUP NDX_FILE=index-wrap.ndx NDX_GROUP=Wrap-AROUND\n")
   out.write("WRAPAROUND ATOMS=wrap-at AROUND=wrap-ar PAIR\n\n")
   # in case of equilibration, add harmonic distance
   # restraints to ordered water molecules
   keys=[]
   for key in ow_dict:
       # define distance CV
       out.write("dr%d: DISTANCE NOPBC ATOMS=%d,%d\n" % (len(keys)+1,key[0],key[1]))
       keys.append(key)
   out.write("\n# Harmonic restraints on ordered waters\n")
   out.write("hr: RESTRAINT STRIDE=%d KAPPA=" % stride)
   for i in range(1,len(keys)):
      out.write("%d," % KAPPA_)
   out.write("%d AT=" % KAPPA_)
   for i in range(0,len(keys)-1):
      out.write("%5.3lf," % ow_dict[keys[i]])
   out.write("%5.3lf ARG=" % ow_dict[keys[-1]])
   for i in range(1,len(keys)):
       out.write("dr%d," % i)
   out.write("dr%d\n" % len(keys))

# add upper walls to all water molecules (ordered + buffer)
# cycle on dictionary
i=1
out.write("\n# CV definition\n")
for key in w_close:
    # define distance CV
    out.write("d%d: DISTANCE NOPBC ATOMS=%d,%d\n" % (i,key,w_close[key]))
    i+=1
# upper limits 
out.write("\n# Upper limits\n")
# define upper wall 
out.write("ul: UPPER_WALLS STRIDE=%d KAPPA=" % stride)
for i in range(1,len(w_close)):
    out.write("%d," % KAPPA_)
out.write("%d AT=" % KAPPA_)
for i in range(1,len(w_close)):
    out.write("%5.3lf," % dmax)
out.write("%5.3lf ARG=" % dmax)
for i in range(1,len(w_close)):
    out.write("d%d," % i)
out.write("d%d\n" % len(w_close))
