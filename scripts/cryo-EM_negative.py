# library to read/write mrc files
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import mrcfile
import numpy as np
import sys
import math
from scipy.stats import pearsonr
from scipy.stats import linregress
from scipy.spatial import distance
from scipy.ndimage import gaussian_filter
import MDAnalysis as mda
from Bio.PDB import *
import argparse

##
## you need the following libraries:
##  conda install -c conda-forge mrcfile mdanalysis biopython
##  conda install -c anaconda scipy
##
## Written by Max Bonomi (mbonomi@pasteur.fr)
##

def get_map_parameters(mrc):
    # initialize dictionary
    mrc_p={}
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
    # return dictionary
    return mrc_p

# auxiliary functions
def indexes2index(i3D, nbin):
    i1D = i3D[2] * nbin[1] * nbin[0] + i3D[1] * nbin[0] + i3D[0]
    return i1D

def index2indexes(i1D, nbin):
    # initialize
    i3D=[0,0,0]
    # calculate
    i3D[0] = int(i1D % nbin[0])
    kk     = int((i1D-i3D[0])/nbin[0])
    i3D[1] = int(kk % nbin[1])
    i3D[2] = int((kk-i3D[1])/nbin[1])
    return tuple(i3D)

# initialize parser
parser = argparse.ArgumentParser(prog='python cryo-EM_negative.py', description='Preprocessing of cryo-EM maps for PLUMED')
parser.add_argument('map', type=str, help='mrc filename')
parser.add_argument('--zone', type=float, default=5.0, metavar='ZONE', help='select voxels close to reference PDB')
parser.add_argument('--zone_PDB', type=str, help='reference PDB for map zoning')
parser.add_argument('--zone_sel', type=str, help='selection for map zoning')

args = parser.parse_args()

#### INPUT
# input MRC file
MRC_=vars(args)["map"]
#### ZONING
do_zoning=False
if(vars(args)["zone_PDB"]!=None):
  do_zoning=True
  ZONE_CUT_ = vars(args)["zone"]
  ZONE_PDB_ = vars(args)["zone_PDB"]
  ZONE_SEL_ = 'all'
  if(vars(args)["zone_sel"]!=None):
    ZONE_SEL_ = vars(args)["zone_sel"]

# open mrc file
mrc = mrcfile.open(MRC_, mode='r+', permissive=True)
# get parameters
mrc_p = get_map_parameters(mrc)
# create a flatten version of the map
data = mrc.data.flatten()
# indexes in the full map (before deletion)
index = np.array(range(0,mrc_p["nbin_tot"]))
print("\n\n General input parameters:")
print("%37s %s" % ("cryo-EM map filename :", MRC_))
print("%37s %d" % ("# voxels :", len(data)))

# if zoning and threshold
if(do_zoning):
  # create MDA universe
  u = mda.Universe(ZONE_PDB_)
  # do atom selection
  at = u.select_atoms(ZONE_SEL_)
  # get positions of selected atoms
  pos = at.positions
  # get minibox around protein
  ilim = []
  for i in range(0,3):
      # min and max coordinate with padding
      xlim = [np.amin(pos[:,i])-ZONE_CUT_, np.amax(pos[:,i])+ZONE_CUT_]
      # and their indexes
      ilim.append([int(math.floor((xlim[0]-mrc_p["x0"][i])/mrc_p["dx"][i])),
                    int(math.ceil((xlim[1]-mrc_p["x0"][i])/mrc_p["dx"][i]))])
  # get mass
  mass = at.total_mass() / 1000.0
  # get resolution
  parser = PDBParser()
  structure = parser.get_structure("SYS", ZONE_PDB_)
  res = structure.header["resolution"]
  # printout
  print("%37s %s" % ("Zoning based on PDB :", ZONE_PDB_))
  print("%37s %s" % ("Zoning atoms selection :", ZONE_SEL_))
  print("%37s %3.1lf" % ("Zoning cutoff [Ang] :", ZONE_CUT_))
  print("%37s %3.2lf" % ("Resolution [Ang] :", res))
  print("%37s %d" % ("Number of atoms :", len(pos)))
  print("%37s %4.1lf" % ("Mass [kDa] :", mass))
  # initialize tokeep list
  tokeep = []
  # cycle on minibox - check boundaries
  for i in range(max(0,ilim[0][0]), min(mrc_p["nbin"][0], ilim[0][1]+1)):
   for j in range(max(0,ilim[1][0]), min(mrc_p["nbin"][1], ilim[1][1]+1)):
    for k in range(max(0,ilim[2][0]), min(mrc_p["nbin"][2], ilim[2][1]+1)):
        # get 1D index
        ii = indexes2index((i,j,k), mrc_p["nbin"])
        # get xyz coord
        xyz=np.array([[mrc_p["x0"][0]+float(i)*mrc_p["dx"][0],
                       mrc_p["x0"][1]+float(j)*mrc_p["dx"][1],
                       mrc_p["x0"][2]+float(k)*mrc_p["dx"][2]]])
        # find minimum distance wrt position array
        mdist = np.amin(distance.cdist(xyz,pos))
        # add to tokeep
        if(mdist<=ZONE_CUT_): tokeep.append(ii)
  # convert tokeep into a numpy array
  tokeep = np.array(tokeep)
else:
  # index of elements to keep
  tokeep = index

# keep only voxels that satisfy the selection criteria (zoning)
data = data[tokeep]
# number of voxels above zero
above = len(data[data>0])
# total number of voxels
nvox = len(data)
# print statistics
print("%37s %4.3lf %4.3lf %4.3lf" % ("Voxel size [Ang] :",mrc_p["dx"][0],mrc_p["dx"][1],mrc_p["dx"][2]))
print("%37s %d" % ("# voxels :", nvox))
print("%37s %d" % ("# voxels above zero :", above))
print("%37s %5.2lf" % ("% voxels above zero :", 100.0*float(above)/float(nvox)))
