# library to read/write mrc files
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import mrcfile
import numpy as np
import sys
import math
import MDAnalysis as mda
import argparse
from scipy.spatial import distance_matrix

def write_to_mrc(bcloud, fname):
    # create empty numpy array
    data = np.zeros((bcloud['NBINS'][0],bcloud['NBINS'][1],bcloud['NBINS'][2]), dtype=np.float32)
    # fill data
    for key in bcloud:
        if type(key) is tuple:
           # add voxels to data array
           data[key[0],key[1],key[2]] = bcloud[key]
    # write to mrc
    with mrcfile.new(fname, overwrite=True) as mrc:
         # set data (we need to transpose the array)
         mrc.set_data(data.T)
         # set voxel size
         mrc.voxel_size = bcloud['DX']
         # set nstart from mrc
         mrc._set_nstart(0,0,0)
         mrc.header.origin.flags.writeable = True
         mrc.header.origin.x = bcloud['XMIN'][0]
         mrc.header.origin.y = bcloud['XMIN'][1]
         mrc.header.origin.z = bcloud['XMIN'][2]
         # update header
         mrc.update_header_from_data()
         mrc.update_header_stats()

def get_map_parameters(mrc):
    # initialize dictionary
    mrc_p={}
    # data organization
    mrc_p["map"] = [mrc.header.mapc, mrc.header.mapr, mrc.header.maps]
    # number of bins
    mrc_p["nbin"] = [mrc.header.nx,mrc.header.ny,mrc.header.nz]
    # total number of bins
    mrc_p["nbin_tot"] = mrc.header.nx*mrc.header.ny*mrc.header.nz
    # origin
    mrc_p["x0"] = [mrc.header.origin.x + float(mrc.header.nxstart) * mrc.voxel_size.x,
                   mrc.header.origin.y + float(mrc.header.nystart) * mrc.voxel_size.y,
                   mrc.header.origin.z + float(mrc.header.nzstart) * mrc.voxel_size.z]
    # dimension of one voxel the in x, y, z directions
    mrc_p["dx"] = [mrc.voxel_size.x,mrc.voxel_size.y,mrc.voxel_size.z]
    # reorder so that is always xyz format
    ijk = [mrc_p["map"].index(1), mrc_p["map"].index(2), mrc_p["map"].index(3)]
    mrc_p["nbin"] = [mrc_p["nbin"][ijk[0]], mrc_p["nbin"][ijk[1]], mrc_p["nbin"][ijk[2]]]
    mrc_p["x0"]   = [mrc_p["x0"][ijk[0]],   mrc_p["x0"][ijk[1]],   mrc_p["x0"][ijk[2]]]
    mrc_p["dx"]   = [mrc_p["dx"][ijk[0]],   mrc_p["dx"][ijk[1]],   mrc_p["dx"][ijk[2]]]
    # return dictionary
    return mrc_p

# initialize parser
parser = argparse.ArgumentParser(prog='python cryo-EM_segment.py', description='Segment a map of a multimer')
parser.add_argument('map', type=str, help='input mrc filename')
parser.add_argument('pdb', type=str, help='pdb filename')
parser.add_argument('sel', type=str, help='pdb selection')
parser.add_argument('--out', type=str, default="output.mrc", help='output mrc filename')

args = parser.parse_args()

#### INPUT
# input MRC file
MRC_=vars(args)["map"]
# PDB file
PDB_=vars(args)["pdb"]
# MDA selection
SEL_=vars(args)["sel"]
# hardcoded slack
SLACK_=5.0

# 1) creating MDA universe
u=mda.Universe(PDB_)
# select atoms for fragmentation
atoms_sel = u.select_atoms("not type H and "+SEL_)
# and their positions
pos_sel = atoms_sel.positions
# select all heavy atoms
atoms_all = u.select_atoms("not type H")
# and their positions
pos_all = atoms_all.positions

# open mrc file
mrc = mrcfile.open(MRC_, mode='r+', permissive=True)
# get parameters
mrc_p = get_map_parameters(mrc)
# get data so that is always zyx format
data = mrc.data.transpose(2-mrc_p["map"].index(3), 2-mrc_p["map"].index(2), 2-mrc_p["map"].index(1))

# get minibox around protein
xmin = np.amin(pos_sel,0)-SLACK_
xmax = np.amax(pos_sel,0)+SLACK_
imin = np.floor((xmin-mrc_p["x0"])/mrc_p["dx"]).astype(int)
imax = np.ceil(( xmax-mrc_p["x0"])/mrc_p["dx"]).astype(int)

# initialize output mrc
output_mrc={'XMIN': mrc_p["x0"],'NBINS': mrc_p["nbin"],'DX': mrc_p["dx"]}

# loop inside minibox
for i in range(max(0,imin[0]), min(mrc_p["nbin"][0], imax[0]+1)):
 for j in range(max(0,imin[1]), min(mrc_p["nbin"][1], imax[1]+1)):
  for k in range(max(0,imin[2]), min(mrc_p["nbin"][2], imax[2]+1)):
      # get xyz coord of voxel
      xyz = [mrc_p["x0"][0]+float(i)*mrc_p["dx"][0],
             mrc_p["x0"][1]+float(j)*mrc_p["dx"][1],
             mrc_p["x0"][2]+float(k)*mrc_p["dx"][2]]
      # calculate distances between xyz and pos_all
      dist = distance_matrix(np.array([xyz]), pos_all)
      # atom closest to the voxel
      iat = np.argmin(dist, axis=1)
      # add voxel to output_mrc if atoms is in the selection
      if(atoms_all[iat] in atoms_sel):
         # add to dict
         output_mrc[(i,j,k)]=data[k][j][i]
 
# write mrc to file
write_to_mrc(output_mrc, vars(args)["out"]) 
