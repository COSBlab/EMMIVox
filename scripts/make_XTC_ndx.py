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
               idx.append(int(riga[i]))
    return idx 

# initialize parser
parser = argparse.ArgumentParser(prog='python make_ndx.py', description='Create GROMACS index.ndx file with atoms selection')
parser.add_argument('group_name',  type=str, nargs='+',  help='Group(s) name(s)')
parser.add_argument('--ndx',       type=str, nargs=1, help='GROMACS ndx file to append output')
args = parser.parse_args()

#### INPUT
# group(s) name(s)
gname = vars(args)["group_name"] 
# output file name
oname = "index.ndx"
if(vars(args)["ndx"]!=None):
  oname=vars(args)["ndx"][0]

# cycle on groups
idx=[]
for g in gname:
    # read indexes from file
    idx += get_atom_indexes(oname,g)
# find maximum index
idx_max = max(idx)

# append to index file
out=open(oname, "a")
out.write("[ System-XTC ]\n")
for i in range(1, idx_max+1):
    out.write("%7d " % i)
    # go to new line
    if(i%15==0 or i==idx_max): out.write("\n")
