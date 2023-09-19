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
import torch

# set device
if(torch.cuda.is_available()):
 device_ = 'cuda'
else:
 device_ = 'cpu'

##
## you need the following libraries:
##  conda install -c conda-forge mrcfile mdanalysis biopython
##  conda install -c anaconda scipy
##
## Written by Max Bonomi (mbonomi@pasteur.fr)
##
## vdw radii [Ang] dictionary
## from ./modules/cctbx_project/cctbx/eltbx/van_der_waals_radii.py
## based on this http://www.phenix-online.org/pipermail/phenixbb/2014-August/021077.html
VDW_={'C': 1.775, 'O': 1.45, 'N': 1.50, 'S': 1.80, 'P': 1.9, 'F': 1.47}
## RPROBE_ and RSHRINK_ from Jiang and Brunger
RPROBE_ = 1.0
RSHRINK_= 1.1
## 5-Gaussians parameters
# map of atom types to A and B coefficients of scattering factor
# f(s) = A * exp(-B*s**2)
# B is in Angstrom squared
# Elastic atomic scattering factors of electrons for neutral atoms
# and s up to 6.0 A^-1: as implemented in PLUMED
### SIGMA
#B_={}
#B_["C"]=np.array([0.1140, 1.0825, 5.4281, 17.8811, 51.1341])
#B_["O"]=np.array([0.0652, 0.6184, 2.9449,  9.6298, 28.2194])
#B_["N"]=np.array([0.0541, 0.5165, 2.8207, 10.6297, 34.3764])
#B_["S"]=np.array([0.0838, 0.7788, 4.3462, 15.5846, 44.63655])
### WEIGHT
#A_={}
#A_["C"]=np.array([0.0489, 0.2091, 0.7537, 1.1420, 0.3555])
#A_["O"]=np.array([0.0365, 0.1729, 0.5805, 0.8814, 0.3121])
#A_["N"]=np.array([0.0267, 0.1328, 0.5301, 1.1020, 0.4215])
#A_["S"]=np.array([0.0915, 0.4312, 1.0847, 2.4671, 1.0852])
# Elastic atomic scattering factors of electrons for neutral atoms
# and s up to 2.0 A^-1: as implemented in phenix
### SIGMA
B_={}
B_["C"]=np.array([0.2465, 1.7100, 6.4094, 18.6113, 50.2523])
B_["O"]=np.array([0.2067, 1.3815, 4.6943, 12.7105, 32.4726])
B_["N"]=np.array([0.2451, 1.7481, 6.1925, 17.3894, 48.1431])
B_["S"]=np.array([0.2681, 1.6711, 7.0267, 19.5377, 50.3888])
B_["P"]=np.array([0.2908, 1.8740, 8.5176, 24.3434, 63.2996])
B_["F"]=np.array([0.2057, 1.3439, 4.2788, 11.3932, 28.7881])
### WEIGHT
A_={}
A_["C"]=np.array([0.0893, 0.2563, 0.7570, 1.0487, 0.3575])
A_["O"]=np.array([0.0974, 0.2921, 0.6910, 0.6990, 0.2039])
A_["N"]=np.array([0.1022, 0.3219, 0.7982, 0.8197, 0.1715])
A_["S"]=np.array([0.2497, 0.5628, 1.3899, 2.1865, 0.7715])
A_["P"]=np.array([0.2548, 0.6106, 1.4541, 2.3204, 0.8477])
A_["F"]=np.array([0.1083, 0.3175, 0.6487, 0.5846, 0.1421])
## create a string of supported atoms for MDAnalysis selection
selt_="type "
for key in B_:
    selt_ += key + " "
## the Gaussian in density (real) space is the FT of scattering factor
## f(r) = A * (pi/B)**1.5 * exp(-pi**2/B*r**2)
## B must be augmented by bfactors (bfact)
## B -> B+bfact/4.0
## we define the following 2 auxiliary variables:
## invs2 = 2.0 * pi**2 / ( B + bfact/4.0 ) 
## pref  = A * (pi/( B + bfact/4.0 ))**1.5
## and we write the forward model as:
## f(r) = pref * exp(-0.5 * invs2 * r**2)
def get_fmod_param(atoms):
    pref=[]; invs2=[]
    # cycle on atoms
    for at in atoms:
        # get atom type
        atype = at.type
        # add to pref list
        pref.append(A_[atype] * pow(math.pi / ( B_[atype] + at.tempfactor/4.0 ), 1.5))
        # add to invs2 list
        invs2.append(2.0 * math.pi * math.pi / ( B_[atype] + at.tempfactor/4.0 ))
    # convert to tensor and copy to device [natoms, 5]
    pref  = torch.tensor(np.array(pref)).to(torch.float).to(device_)
    invs2 = torch.tensor(np.array(invs2)).to(torch.float).to(device_)
    return pref, invs2

def get_map_parameters(mrc):
    # initialize dictionary
    mrc_p={}
    # data organization
    mrc_p["map"] = [mrc.header.mapc, mrc.header.mapr, mrc.header.maps]
    # number of bins
    mrc_p["nbin"] = [int(mrc.header.nx),int(mrc.header.ny),int(mrc.header.nz)]
    # total number of bins
    mrc_p["nbin_tot"] = int(mrc.header.nx*mrc.header.ny*mrc.header.nz)
    # origin
    mrc_p["x0"] = [float(mrc.header.origin.x + float(mrc.header.nxstart) * mrc.voxel_size.x),
                   float(mrc.header.origin.y + float(mrc.header.nystart) * mrc.voxel_size.y),
                   float(mrc.header.origin.z + float(mrc.header.nzstart) * mrc.voxel_size.z)]
    # dimension of one voxel the in x, y, z directions
    mrc_p["dx"] = [float(mrc.voxel_size.x),float(mrc.voxel_size.y),float(mrc.voxel_size.z)]
    # reorder so that is always xyz format
    ijk = [mrc_p["map"].index(1), mrc_p["map"].index(2), mrc_p["map"].index(3)]
    mrc_p["nbin"] = [mrc_p["nbin"][ijk[0]], mrc_p["nbin"][ijk[1]], mrc_p["nbin"][ijk[2]]]
    mrc_p["x0"]   = [mrc_p["x0"][ijk[0]],   mrc_p["x0"][ijk[1]],   mrc_p["x0"][ijk[2]]]
    mrc_p["dx"]   = [mrc_p["dx"][ijk[0]],   mrc_p["dx"][ijk[1]],   mrc_p["dx"][ijk[2]]]
    # return dictionary
    return mrc_p

def index2indexes(i1D, nbin):
    # initialize
    i3D=[0,0,0]
    # calculate
    i3D[0] = int(i1D % nbin[0])
    kk     = int((i1D-i3D[0])/nbin[0])
    i3D[1] = int(kk % nbin[1])
    i3D[2] = int((kk-i3D[1])/nbin[1])
    return tuple(i3D)

# get vdw radii
def get_vdw_radii(atoms):
    vdw=[]
    for at in atoms:
        vdw.append(VDW_[at.type])
    return np.array(vdw)

# select voxels around PDB
def get_voxels_from_PDB(at,mrc_p,data,thres,rprobe,rshrink):
  # get non-water atoms
  atoms = at.select_atoms("not resname HOH TIP TIP3 SOL")
  # positions tensor on device [nat, 3]
  pos_g = torch.tensor(atoms.positions).to(torch.float).to(device_)
  # get vdw radii
  vdw = get_vdw_radii(atoms)
  # vdw radii: convert to tensor and copy to device [nat]
  vdw_g = torch.tensor(vdw).to(torch.float).to(device_)
  # put map on GPU [nz, ny, nz]
  data_g = torch.tensor(data).to(torch.float).to(device_)
  # put map info on GPU
  x0   = torch.tensor(mrc_p["x0"]).to(torch.float).to(device_)
  dx   = torch.tensor(mrc_p["dx"]).to(torch.float).to(device_)
  nbin = torch.tensor(mrc_p["nbin"]).to(torch.int).to(device_)
  # get minibox around atoms
  xmin = torch.amin(pos_g, 0)-rprobe-torch.amax(vdw_g, 0)
  xmax = torch.amax(pos_g, 0)+rprobe+torch.amax(vdw_g, 0)
  # indices for slicing full matrix
  imin = torch.max(torch.floor((xmin-x0)/dx).int(), torch.tensor([0,0,0]).to(torch.int).to(device_))
  imax = torch.min(torch.ceil( (xmax-x0)/dx).int()+1, nbin)
  # indexes of entries above threshold in minibox (tuple of 3 tensors)
  ind = torch.where(data_g[imin[2]:imax[2],imin[1]:imax[1],imin[0]:imax[0]]>thres)
  # indexes in full/original map (tuple of 3 tensors)
  ind3D_l = (imin[0]+ind[2], imin[1]+ind[1], imin[2]+ind[0])
  # 1D indexes in full/original (flattened) map
  iivox = nbin[1] * nbin[0] * ind3D_l[2] + nbin[0] * ind3D_l[1] + ind3D_l[0]
  # 3D indexes in full/original map [nvox, 3]
  ind3D = torch.cat((ind3D_l[0][:,None], ind3D_l[1][:,None], ind3D_l[2][:,None]), dim=1)
  # voxel positions [nvox, 3]
  vox_g = x0 + dx * ind3D.float()
  # total number of voxels
  nvox = len(vox_g)
  # divide in chunk
  nchunk = 10000
  niter = max(1, int(math.ceil(nvox / nchunk)))
  # initialize indexes tensor [nvox]
  iivox_g = torch.arange(start=0, end=nvox, dtype=torch.long, device=device_)
  # cycle over chunks
  for i in range(0, niter):
      # boundary
      i0 = i * nchunk
      i1 = i0 + nchunk
      # last iteraction
      if(i==niter-1): i1 = nvox
      # check if enough data
      if(i1-i0<=1): continue
      # calculate distances between atom positions and chunk of voxels
      dist = torch.cdist(pos_g, vox_g[i0:i1,:])-vdw_g[:,None].expand(-1, i1-i0)
      # find voxels with neighboring atoms
      cut = torch.any(torch.lt(dist, rprobe), 0)
      # create or add to list of voxels to retain
      if(i==0):
         tokeep_g = iivox_g[i0:i1][cut]
      else:
         tokeep_g = torch.cat((tokeep_g, iivox_g[i0:i1][cut]), 0)
  # create a tensor with discarded voxels
  mask_g = torch.ones(nvox, dtype=torch.bool, device=device_)
  # these are not solvent voxels
  mask_g[tokeep_g] = False
  # indices of voxels to remove (solvent mask)
  toremove_g = iivox_g[mask_g]
  # total number of voxels in solvent mask
  nvox = len(toremove_g)
  # divide in small chunks
  nchunk = 100
  niter = max(1, int(math.ceil(nvox / nchunk)))
  # cycle over chunks
  for i in range(0, niter):
      # boundary
      i0 = i * nchunk
      i1 = i0 + nchunk
      # last iteration
      if(i==niter-1): i1 = nvox
      # check if enough data
      if(i1-i0<=1): continue
      # calculate distances between voxels to keep and chunk of solvent mask
      dist = torch.cdist(vox_g[tokeep_g], vox_g[toremove_g[i0:i1]])
      # find voxels within shrink radius from solvent mask
      if(i==0):
         cut = torch.all(torch.ge(dist, rshrink), 1)
      else:
         cut = torch.logical_and(cut, torch.all(torch.ge(dist, rshrink), 1))
  # tensor of indices in full/original map + move to CPU as list
  tokeep = iivox[tokeep_g[cut]].cpu().tolist()
  return tokeep

def get_voxels(tookeep, mrc_p):
    vox=[]
    for ii in tokeep:
        # get coordinates
        ijk = index2indexes(ii,mrc_p["nbin"])
        ## get voxel position
        xyz=[mrc_p["x0"][0]+float(ijk[0])*mrc_p["dx"][0],
             mrc_p["x0"][1]+float(ijk[1])*mrc_p["dx"][1],
             mrc_p["x0"][2]+float(ijk[2])*mrc_p["dx"][2]]
        # add to list
        vox.append(xyz)
    # convert to tensor and copy to device [nvox, 3]
    vox_g = torch.tensor(vox).to(torch.float).to(device_)
    return vox_g

def get_model_map(atoms,pref,invs2,vox_g,nvox,cut):
    # convert to tensor and copy to device [nat, 3]
    pos_g = torch.tensor(atoms.positions).to(torch.float).to(device_)
    # total number of voxels
    nvox = len(vox_g)
    # divide in chunk
    nchunk = 10000
    niter = max(1, int(math.ceil(nvox / nchunk)))
    # cycle over chunks
    for i in range(niter):
        # boundary
        i0 = i * nchunk
        i1 = i0 + nchunk
        # last iteraction
        if(i==niter-1): i1 = nvox
        # check if enough data
        if(i1-i0<=1): continue
        # create or add to tensor of close atoms (distance lower cutoff)
        if(i==0):
           nl_ia, nl_iv = torch.where(torch.cdist(pos_g, vox_g[i0:i1]) < cut)
        else:
           nl_ia_t, nl_iv_t = torch.where(torch.cdist(pos_g, vox_g[i0:i1]) < cut)
           # cat to full tensor
           nl_ia = torch.cat((nl_ia, nl_ia_t), 0)
           nl_iv = torch.cat((nl_iv, nl_iv_t+i0), 0)
    # get positions atoms in neighbor list [nvox_nl, 3]
    pos_nl = torch.index_select(pos_g,0,nl_ia)
    # get voxels in neighbor list [nvox_nl, 3]
    vox_nl = torch.index_select(vox_g,0,nl_iv)
    # calculate distance squared positions / voxels in neighbor list [nvox_nl]
    dist2 = torch.sum(torch.pow(pos_nl-vox_nl,2), 1)
    # expand distances cutoff [nvox_nl, 5]
    dist2 = dist2[:,None].expand(-1,5)
    # calculate map [nvox_nl]
    map_g = torch.sum(torch.index_select(pref,0,nl_ia) * torch.exp(-0.5 * torch.index_select(invs2,0,nl_ia) * dist2), 1)
    # initialize sum map on device [nvox]
    map_sum = torch.zeros(nvox, dtype=torch.float, device=device_)
    # sum entries of map_g [nvox]
    map_sum.index_add_(0, nl_iv, map_g)
    # and return map [on GPU]
    return map_sum

# initialize parser
parser = argparse.ArgumentParser(prog='python cryo-EM_validate.py', description='Validation of structural models and trajectories against cryo-EM map')
parser.add_argument('map', type=str, help='mrc filename')
parser.add_argument('--threshold', type=float, metavar='TH', default=-10000.0, help='map thresold')
parser.add_argument('--pdbA', type=str, help='name of PDB A')
parser.add_argument('--selA', type=str, help='selection for PDB A')
parser.add_argument('--pdbB', type=str, help='name of PDB B')
parser.add_argument('--selB', type=str, help='selection for PDB B')
parser.add_argument('--pdbC', type=str, help='name of PDB C')
parser.add_argument('--trjC', type=str, help='name of trajectory C')
parser.add_argument('--selC', type=str, help='selection for PDB C')
parser.add_argument('--cut',  type=float, default=4.0, metavar='CUT', help='cutoff Gaussians')

args = parser.parse_args()

#### INPUT
# input MRC file
MRC_=vars(args)["map"]
# thresold
THRES_=vars(args)["threshold"]
#### PDB A
if(vars(args)["pdbA"]==None):
  print("Please specify at least one PDB name with --pdbA")
  exit()
PDBA_ = vars(args)["pdbA"]
# and selection for PDB A
SELA_ = selt_
if(vars(args)["selA"]!=None):
   SELA_ += "and (" + vars(args)["selA"] + ")"
#### PDB B
do_pdbB = False
if(vars(args)["pdbB"]!=None):
   do_pdbB = True
   PDBB_ = vars(args)["pdbB"]
   SELB_ = selt_
   if(vars(args)["selB"]!=None):
      SELB_ += "and (" + vars(args)["selB"] + ")"
#### PDB C
do_pdbC = False
if(vars(args)["pdbC"]!=None):
   do_pdbC = True
   PDBC_ = vars(args)["pdbC"]
   SELC_ = selt_
   if(vars(args)["selC"]!=None):
      SELC_ += "and (" + vars(args)["selC"] + ")"
   if(vars(args)["trjC"]==None):
      print("Please specify a trajectory file with --trjC")
      exit()
   TRJC_ = vars(args)["trjC"]

#### GAUSSIAN CUTOFF
CUT_ = vars(args)["cut"]

# open mrc file
mrc = mrcfile.open(MRC_, mode='r+', permissive=True)
# get parameters
mrc_p = get_map_parameters(mrc)
# get data so that is always zyx format
data = mrc.data.transpose(2-mrc_p["map"].index(3), 2-mrc_p["map"].index(2), 2-mrc_p["map"].index(1))

# prepare voxels list
tokeepA=[]; tokeepB=[]; tokeepC=[]

### PDB A
# create MDA universe
uA = mda.Universe(PDBA_)
# do atom selection
atA = uA.select_atoms(SELA_)
# get list of voxels to keep
tokeepA = get_voxels_from_PDB(atA,mrc_p,data,THRES_,RPROBE_,RSHRINK_)

### PDB B
if(do_pdbB):
  # create MDA universe
  uB = mda.Universe(PDBB_)
  # do atom selection
  atB = uB.select_atoms(SELB_)
  # get list of voxels to keep
  tokeepB = get_voxels_from_PDB(atB,mrc_p,data,THRES_,RPROBE_,RSHRINK_)

### PDB C
if(do_pdbC):
  # create MDA universe
  uC = mda.Universe(PDBC_,TRJC_)
  nframes = len(uC.trajectory)
  # do atom selection
  atC = uC.select_atoms(SELC_)
  # cycle on frames
  #for i in range(0, nframes, 50):
  #    # go to frame i
  #    uC.trajectory[i]
  #    # add to list of voxels to keep
  #    tokeepC += get_voxels_from_PDB(atC,mrc_p,data,THRES_,RPROBE_,RSHRINK_)
  # remove duplicates
  #tokeepC = list(set(tokeepC))

# get final list of voxels
tokeep = list(set(tokeepA+tokeepB+tokeepC))
# get voxel positions on GPU
vox_g = get_voxels(tokeep, mrc_p)
# number of voxels
nvox = len(tokeep)
# create a flatten version of the map
data = data.flatten()

## PRINT OUT SOME INFO
print("\n\n General input parameters:")
print("   * Model A:")
print("%37s %s" % ("PDB filename :", PDBA_))
print("%37s %s" % ("Selection :", SELA_))
print("%37s %d" % ("# atoms :", len(atA)))
if(do_pdbB):
  print("   * Model B:")
  print("%37s %s" % ("PDB filename :", PDBB_))
  print("%37s %s" % ("Selection :", SELB_))
  print("%37s %s" % ("# atoms :", len(atB)))
if(do_pdbC):
  print("   * Model C:")
  print("%37s %s" % ("PDB filename :", PDBC_))
  print("%37s %s" % ("Selection :", SELC_))
  print("%37s %s" % ("# atoms :", len(atC)))
  print("%37s %s" % ("Trajectory filename :", TRJC_))
  print("%37s %d" % ("Trajectory # frames :", nframes))
print("   * Cryo-EM map:")
print("%37s %s" % ("Filename :", MRC_))
print("%37s %d" % ("# voxels :", len(data)))
print("%37s %lf" % ("Threshold :", THRES_))
print("%37s %d" % ("# voxels for cc calculation:", nvox))
print("%37s %lf" % ("Gaussian cutoff :", CUT_))
print("\n")
print("Cross-correlations:")

# keep only voxels that satisfy the selection criteria
data = data[np.array(tokeep)]

## MODEL A
# prepare fmod parameters
prefA, invs2A = get_fmod_param(atA)
# get model map A [still on GPU]
mapA = get_model_map(atA,prefA,invs2A,vox_g,nvox,CUT_)
# do linear regression mapA vs data
slope, intercept, r_value, p_value, std_err = linregress(data, mapA.cpu().numpy())
# print correlation
print("%37s %lf" % ("CC_mask+ mapA / data :", r_value))

## MODEL B
if(do_pdbB):
   # prepare fmod parameters
   prefB, invs2B = get_fmod_param(atB)
   # get model map B [still on GPU]
   mapB = get_model_map(atB,prefB,invs2B,vox_g,nvox,CUT_)
   # do linear regression mapB vs data
   slope, intercept, r_value, p_value, std_err = linregress(data, mapB.cpu().numpy())
   # print correlation
   print("%37s %lf" % ("CC_mask+ mapB / data :", r_value))

## MODEL C
if(do_pdbC):
   # prepare fmod parameters
   prefC, invs2C = get_fmod_param(atC)
   # prepare map on GPU
   mapC = torch.zeros(nvox, dtype=torch.float, device=device_)
   # cycle on frames
   for i in range(0, nframes):
       # load frame
       uC.trajectory[i]
       # add to mapC on GPU
       mapC += get_model_map(atC,prefC,invs2C,vox_g,nvox,CUT_)
   # calculate average map
   mapC /= float(nframes)
   # do linear regression mapC vs data
   slope, intercept, r_value, p_value, std_err = linregress(data, mapC.cpu().numpy())
   # print correlation
   print("%37s %lf" % ("CC_mask+ mapC / data :", r_value))
