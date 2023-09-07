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

# initialize parser
parser = argparse.ArgumentParser(prog='python error-residue.py', description='Get map error per residue')
parser.add_argument('in_pdb',  type=str, help='input PDB file')
parser.add_argument('map',     type=str, help='input map in PLUMED format')
parser.add_argument('out',     type=str, help='output error file')
parser.add_argument('--sel',   type=str, nargs=1, help="MDAnalysis selection for input PDB")
args = parser.parse_args()

#### INPUT
# input PDB
IN_PDB=vars(args)["in_pdb"]
# map in PLUMED format
MAP_=vars(args)["map"]
# output error file
OUT_=vars(args)["out"]
# selection
SEL="not type H"
if(vars(args)["sel"]!=None):
  SEL=vars(args)["sel"][0]
# cutoff
CUT_=4.0

# create universe
pdb=mda.Universe(IN_PDB)
# select atoms
sel=pdb.select_atoms(SEL)
# separate selection for backbone
sel_back = sel.select_atoms('backbone')

# for each atom, store: position, resid, chainid
pos=[]; resid=[]; segid=[]
for at in sel:
    pos.append(at.position)
    resid.append(at.resid)
    segid.append(at.segid)
# convert position list to numpy array
pos = np.array(pos)

# prepare error dictionary:
# key: pair of residue,chain ids
# values: list of relative errors of the voxels associated to key
data_err={}; data_d_back={}; data_d_side={}
# cycle on voxels
for lines in open(MAP_, "r").readlines():
    riga=lines.strip().split()
    # skip comment line
    if(riga[0]=="#!"): continue
    # get voxel center - in Angstrom
    pos_v = 10.0*np.array((float(riga[1]),float(riga[2]),float(riga[3])))
    # get density
    d = float(riga[4])
    # get relative error
    err = float(riga[5]) / d
    # calculate distances between voxel center and all atoms
    all_dist = np.linalg.norm(np.array([pos_v])-pos, axis=1)
    # index of closest atom
    i_min = all_dist.argmin()
    # if within cutoff
    if(all_dist[i_min]<CUT_):
       # define key in data_err
       key=(resid[i_min],segid[i_min])
       # add to dictionary
       if(key in data_err):
         data_err[key].append(err)
       else:
         data_err[key]=[err]
       # add to density dictionaries based on sidechain/backbone
       # this is a backbone atom
       if(sel[i_min] in sel_back):
            if(key in data_d_back):
               data_d_back[key].append(d)
            else:
               data_d_back[key]=[d]
       # this is a sidechain atom
       else:
            if(key in data_d_side):
               data_d_side[key].append(d)
            else:
               data_d_side[key]=[d]

# prepare output file
out=open(OUT_, "w")
# write header
out.write("#! Rs  Ch    MedErr   MedDBack  MedDSid\n")
# cycle on chains
for ch in sel.segments:
  # cycle on residues
  for r in ch.residues:
      # dictionary key
      key=(r.resid,r.segid)
      # get array of errors
      err = np.array([0.0])
      # convert to numpy array
      if key in data_err:
         err = np.array(data_err[key])
      # find median
      err_med = np.median(err)
      # get array of backbone density 
      d_back = np.array([0.0])
      # convert to numpy array
      if key in data_d_back:
         d_back = np.array(data_d_back[key])
      # find median
      d_back_med = np.median(d_back)
      # get array of sidechain density 
      d_side = np.array([0.0])
      # convert to numpy array
      if key in data_d_side:
         d_side = np.array(data_d_side[key])
      # find median
      d_side_med = np.median(d_side)
      # print to file
      out.write("%5d %3s  %8.3lf   %8.3lf %8.3lf\n" % (key[0],key[1],err_med,d_back_med,d_side_med))
