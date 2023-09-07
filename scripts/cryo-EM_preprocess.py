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

# write numpy 3D array to mrc
def write_data_to_mrc(filename, data, mrc):
    # creating new mrc file
    mrc_map = mrcfile.new(filename, overwrite=True)
    # set headers from mrc
    mrc_map._set_nstart(mrc.header.nxstart,mrc.header.nystart,mrc.header.nzstart)
    mrc_map.header.origin.flags.writeable = True
    mrc_map.header.origin.x = mrc.header.origin.x
    mrc_map.header.origin.y = mrc.header.origin.y
    mrc_map.header.origin.z = mrc.header.origin.z
    mrc_map.header.cella.flags.writeable = True
    mrc_map.header.cella.x = mrc.header.cella.x
    mrc_map.header.cella.y = mrc.header.cella.y
    mrc_map.header.cella.z = mrc.header.cella.z
    # set data
    mrc_map.set_data(data)
    # update header from data
    mrc_map.update_header_from_data()
    # close file
    mrc_map.close()

# initialize parser
parser = argparse.ArgumentParser(prog='python cryo-EM_preprocess.py', description='Preprocessing of cryo-EM maps for PLUMED')
parser.add_argument('map', type=str, help='mrc filename')
parser.add_argument('cutoff', type=float, help='correlation cutoff')
parser.add_argument('plumed', type=str, help='PLUMED output map file name')
parser.add_argument('--threshold', type=float, metavar='TH', default=0.0, help='map thresold')
parser.add_argument('--box', type=float, default=2.0, help='side of minibox to calculate correlation')
parser.add_argument('--halfmaps', type=str, nargs=2, metavar='HMAP', help='mrc filenames of the two half-maps')
parser.add_argument('--smooth', type=float, default=1.0, metavar='NVOX', help='apply Gaussian filtering (smoothing) to error map')
parser.add_argument('--errmap', type=str, help='output error map mrc file')
parser.add_argument('--zone', type=float, default=3.5, metavar='ZONE', help='select voxels close to reference PDB')
parser.add_argument('--zone_PDB', type=str, help='reference PDB for map zoning')
parser.add_argument('--zone_sel', type=str, help='selection for map zoning')
parser.add_argument('--double', default=False, action='store_true')

args = parser.parse_args()

#### INPUT
# input MRC file
MRC_=vars(args)["map"]
# thresold
THRES_=vars(args)["threshold"]
# half side of cubic box to calculate correlation in Angstrom
# typically 2 Ang
CLMAX_=vars(args)["box"]
# pearson cutoff
# typically 0.8
PCUT_=vars(args)["cutoff"]
#### OUTPUT
# plumed file
PLUMED_=vars(args)["plumed"]
#### HALF MAPS / ERROR MAP
do_halfmaps=False; do_smooth=False; do_errmap=False
if(vars(args)["halfmaps"]!=None):
  do_halfmaps=True
  MRC1_ = vars(args)["halfmaps"][0]
  MRC2_ = vars(args)["halfmaps"][1]
  SMOOTH_ = vars(args)["smooth"]
  if(SMOOTH_>0.): do_smooth=True
  ERRMAP_ = vars(args)["errmap"]
  if(ERRMAP_!=None): do_errmap=True
#### ZONING
do_zoning=False
if(vars(args)["zone_PDB"]!=None):
  do_zoning=True
  ZONE_CUT_ = vars(args)["zone"]
  ZONE_PDB_ = vars(args)["zone_PDB"]
  ZONE_SEL_ = 'all'
  if(vars(args)["zone_sel"]!=None):
    ZONE_SEL_ = vars(args)["zone_sel"]
### use double precision
if(vars(args)["double"]):
 torch.set_default_dtype(torch.double)
else:
 torch.set_default_dtype(torch.float)

# open mrc file
mrc = mrcfile.open(MRC_, mode='r+', permissive=True)
# get parameters
mrc_p = get_map_parameters(mrc)
# get data so that is always zyx format
data_c = mrc.data.transpose(2-mrc_p["map"].index(3), 2-mrc_p["map"].index(2), 2-mrc_p["map"].index(1))
# put map on GPU [nz, ny, nz]
data_g = torch.tensor(data_c).to(device_)
# put map info on GPU
x0   = torch.tensor(mrc_p["x0"]).to(device_)
dx   = torch.tensor(mrc_p["dx"]).to(device_)
nbin = torch.tensor(mrc_p["nbin"]).to(torch.int).to(device_)

# printout
print("\n\n General input parameters:")
print("%37s %s" % ("cryo-EM map filename :", MRC_))
print("%37s %d" % ("# voxels :", len(data_c.flatten())))

# if zoning and threshold
if(do_zoning):
  # create MDA universe
  u = mda.Universe(ZONE_PDB_)
  # do atom selection
  atoms = u.select_atoms(ZONE_SEL_)
  # positions tensor on device [nat, 3]
  pos_g = torch.tensor(atoms.positions).to(x0.dtype).to(device_)
  # get mass
  mass = atoms.total_mass() / 1000.0
  # get minibox around atoms
  xmin = torch.amin(pos_g, 0)-ZONE_CUT_
  xmax = torch.amax(pos_g, 0)+ZONE_CUT_
  # indices for slicing full matrix
  imin = torch.max(torch.floor((xmin-x0)/dx).int(), torch.tensor([0,0,0]).to(torch.int).to(device_))
  imax = torch.min(torch.ceil( (xmax-x0)/dx).int()+1, nbin)
  # get resolution
  parser = PDBParser()
  structure = parser.get_structure("SYS", ZONE_PDB_)
  res = structure.header["resolution"]
  # printout
  print("%37s %s" % ("Zoning based on PDB :", ZONE_PDB_))
  print("%37s %s" % ("Zoning atoms selection :", ZONE_SEL_))
  print("%37s %3.1lf" % ("Zoning cutoff [Ang] :", ZONE_CUT_))
  print("%37s %3.2lf" % ("Resolution [Ang] :", res))
  print("%37s %d" % ("Number of atoms :", len(pos_g)))
  print("%37s %4.1lf" % ("Mass [kDa] :", mass))
  # indexes of entries above threshold in minibox (tuple of 3 tensors)
  ind = torch.where(data_g[imin[2]:imax[2],imin[1]:imax[1],imin[0]:imax[0]]>THRES_)
  # indexes in full/original map (tuple of 3 tensors)
  ind3D_l = (imin[0]+ind[2], imin[1]+ind[1], imin[2]+ind[0])
  # 1D indexes in full/original (flattened) map
  iivox = nbin[1] * nbin[0] * ind3D_l[2] + nbin[0] * ind3D_l[1] + ind3D_l[0]
  # 3D indexes in full/original map [nvox, 3]
  ind3D = torch.cat((ind3D_l[0][:,None], ind3D_l[1][:,None], ind3D_l[2][:,None]), dim=1)
  # voxel positions [nvox, 3]
  vox_g = x0 + dx * ind3D
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
      # calculate distances between atom positions and chunk of voxels
      cut = torch.any(torch.le(torch.cdist(vox_g[i0:i1,:], pos_g), ZONE_CUT_), 1)
      # create or add to list of voxels to retain
      if(i==0):
         tokeep_g = iivox[i0:i1][cut]
      else:
         tokeep_g = torch.cat((tokeep_g, iivox[i0:i1][cut]), 0)
else:
  # indexes of entries above threshold (tuple of 3 tensors)
  ind = torch.where(data_g>THRES_)
  # indexes in full/original map (tuple of 3 tensors)
  ind3D_l = (ind[2], ind[1], ind[0])
  # 1D indexes in full/original (flattened) map
  tokeep_g = nbin[1] * nbin[0] * ind3D_l[2] + nbin[0] * ind3D_l[1] + ind3D_l[0]

# convert tokeep tensor to numpy array
index = tokeep_g.cpu().numpy()
# keep only voxels that satisfy the selection criteria (threshold/zoning)
data = data_c.flatten()[index]
# index for sorting data
index_srt = np.argsort(data)[::-1]
# now sort both data and index
data  = data[index_srt]
index = index[index_srt]

# total number of voxels
nvox = len(data)
# calculate integral of density on filtered map
normd = np.sum(data)*mrc_p["dx"][0]*mrc_p["dx"][1]*mrc_p["dx"][2]
# calculate median density
med = np.median(data)
# number of voxels above/below median
above_med = len(data[data>med])
below_med = nvox-above_med
# print statistics
print("%37s %4.3lf %4.3lf %4.3lf" % ("Voxel size [Ang] :",mrc_p["dx"][0],mrc_p["dx"][1],mrc_p["dx"][2]))
print("%37s %lf" % ("Threshold :", THRES_))
print("%37s %d"  % ("# filtered voxels :", nvox))
print("%37s %lf" % ("Integral of density (NORM_DENSITY):", normd))
print("%37s %lf" % ("Median density :", med))
print("%37s %d / %d" % ("# voxels above/below median :", above_med, below_med))
if(do_zoning):
  print("%37s %lf" % ("# voxels per kDa :", float(nvox)/mass))
  print("%37s %lf" % ("# voxels per atom :", float(nvox)/float(len(pos_g))))
# prepare list of voxels for PLUMED:
# kcenters contains the position in the 3D grid
# dcenters contains the value of the density
kcenters=[]; dcenters=[];

# determining box dimensions to calculate local correlation
NMAX_=[]
for i in range(3):
 NMAX_.append(int(math.ceil(CLMAX_/mrc_p["dx"][i])))

print("\n Creating correlation-filtered map:")
print("%37s %4.3lf" % ("  Minibox size [Ang] :", CLMAX_))
print("%37s %4.3lf" % ("  Correlation cutoff :", PCUT_))

# in case of no filtering
if(PCUT_>=1.0):
  for i in range(len(data)):
      # calculate index 3D
      ii = index2indexes(index[i], mrc_p["nbin"])
      # add to list of voxels
      dcenters.append(data[i]) # density
      kcenters.append(ii)      # index 3D
# otherwise...
else:
  # prepare boolean list to skip voxels
  toskip = len(data_c.flatten()) * [False]
  # loop over voxels
  for i in range(len(data)):
     # 0) skip voxel if marked as correlated
     if toskip[index[i]]: continue
     # 1) value of the density map
     dd = data[i]
     # get 3D index
     ii = index2indexes(index[i],mrc_p["nbin"])
     # add to list of voxels
     dcenters.append(dd) # density
     kcenters.append(ii) # index 3D
     # define indexes for minibox used to evaluate cc
     i0, i1 = max(0, ii[0]-NMAX_[0]), min(mrc_p["nbin"][0], ii[0]+NMAX_[0]+1)
     j0, j1 = max(0, ii[1]-NMAX_[1]), min(mrc_p["nbin"][1], ii[1]+NMAX_[1]+1)
     k0, k1 = max(0, ii[2]-NMAX_[2]), min(mrc_p["nbin"][2], ii[2]+NMAX_[2]+1)
     # calculate max displacement (needed at borders)
     nmax0 = min(ii[0]-i0, i1-ii[0]-1)
     nmax1 = min(ii[1]-j0, j1-ii[1]-1)
     nmax2 = min(ii[2]-k0, k1-ii[2]-1)
     # 2) now calculate correlation coefficient in a cube of side 2*NMAX_+1 centered on ii
     # loop over displacements
     for c0 in range(nmax0+1):
      for c1 in range(nmax1+1):
       for c2 in range(nmax2+1):
            # get slices
            data1 = data_c[k0:k1-c2, j0:j1-c1, i0:i1-c0].flatten()
            data2 = data_c[k0+c2:k1, j0+c1:j1, i0+c0:i1].flatten()
            # pearson correlation coefficient
            pc = pearsonr(data1, data2)[0]
            # remove voxels
            if(pc>PCUT_):
               for di in set([-c0,c0]):
                for dj in set([-c1,c1]):
                 for dk in set([-c2,c2]):
                    # get 3D indices
                    i3D = (ii[0]+di, ii[1]+dj, ii[2]+dk)
                    # get 1D index
                    i1D = indexes2index(i3D, mrc_p["nbin"])
                    # mark toskip
                    toskip[i1D] = True

# total number of voxels after cc-filtering
nvoxr = len(dcenters)
# print log
print("%37s %d" % ("  Reduced number of voxels :", nvoxr))
# statistics on reduced map
datar=np.array(dcenters)
# number of voxels above/below (old) median
above_medr = len(datar[datar>med])
below_medr = nvoxr-above_medr
print("%37s %d / %d" % ("  # voxels above/below median :", above_medr, below_medr))
# discarding percentages
p = float(nvox - nvoxr) / float(nvox) * 100.0
p_a = 0.0; p_b = 0.0
if(nvox > nvoxr):
   p_a = float( above_med - above_medr ) / float(nvox - nvoxr) * 100.0
   p_b = float( below_med - below_medr ) / float(nvox - nvoxr) * 100.0
print("%37s %4.1lf" % ("  % discard :", p))
print("%37s %4.1lf / %4.1lf" % ("  % discard above/below median :", p_a, p_b))
if(do_zoning):
  print("%37s %lf" % ("# voxels per kDa :", float(nvoxr)/mass))
  print("%37s %lf" % ("# voxels per atom :", float(nvoxr)/float(len(pos_g))))
#
# do half map analysis
if(do_halfmaps):
 print("\n Doing half-maps analysis:")
 print("%37s %s" % ("half-map1 filename :",MRC1_))
 print("%37s %s" % ("half-map2 filename :",MRC2_))
 # open mrc file 1
 mrc1 = mrcfile.open(MRC1_, permissive=True)
 # get parameters
 mrc1_p = get_map_parameters(mrc1)
 # open mrc file 2
 mrc2 = mrcfile.open(MRC2_, permissive=True)
 # get parameters
 mrc2_p = get_map_parameters(mrc2)
 # check parameters: comparison with full map
 for key in mrc_p:
     if(mrc_p[key]!=mrc1_p[key] or mrc_p[key]!=mrc2_p[key]):
        print(" ERROR: half maps should be on the same grid as full map!")
        exit()
 # get data so that is always zyx format
 data1_c = mrc1.data.transpose(2-mrc_p["map"].index(3), 2-mrc_p["map"].index(2), 2-mrc_p["map"].index(1))
 data2_c = mrc2.data.transpose(2-mrc_p["map"].index(3), 2-mrc_p["map"].index(2), 2-mrc_p["map"].index(1))
 # remove filtered voxels before linear regression
 data  = data_c.flatten()[index]
 data1 = data1_c.flatten()[index]
 data2 = data2_c.flatten()[index]
 # do linear regression with full map
 slope1, intercept1, r_value1, p_value1, std_err1 = linregress(data1,data)
 slope2, intercept2, r_value2, p_value2, std_err2 = linregress(data2,data)
 # calculate correlation (Chimera-like)
 cc1 = np.dot(data, slope1*data1+intercept1) / np.linalg.norm(data) / np.linalg.norm(slope1*data1+intercept1)
 cc2 = np.dot(data, slope2*data2+intercept2) / np.linalg.norm(data) / np.linalg.norm(slope2*data2+intercept2)
 print("%37s %4.3lf" % ("Correlation half-map1 to full map :",cc1))
 print("%37s %4.3lf" % ("Correlation half-map2 to full map :",cc2))
 # error map
 err = 0.5 * np.sqrt(np.power(slope1*data1_c+intercept1-data_c,2)+np.power(slope2*data2_c+intercept2-data_c,2))
 # smoothing with neighbors
 if(do_smooth):
    print("\n Applying Gaussian filtering (smoothing) to error map:")
    print("%37s %2.1lf" % ("sigma [# voxels] :",SMOOTH_))
    err = gaussian_filter(err, sigma=SMOOTH_*mrc_p["dx"][0])
 # set error to zero in deleted voxels
 index_tmp = np.array(range(mrc_p["nbin_tot"]))
 todelete = np.delete(index_tmp, index)
 for v in todelete:
     ii = index2indexes(v, mrc_p["nbin"])
     err[ii[2],ii[1],ii[0]] = 0.0
 # printing error map to mrc file
 if(do_errmap):
    print("\n Writing error map to mrc file:")
    print("%37s %s" % ("error map filename :",ERRMAP_))
    write_data_to_mrc(ERRMAP_, err, mrc)
    # calculate relative error map
    err_rel = np.copy(err)
    # cycle on voxels to keep
    for v in index:
        ii = index2indexes(v, mrc_p["nbin"])
        err_rel[ii[2],ii[1],ii[0]] /= data_c[ii[2],ii[1],ii[0]]
    # write to file
    write_data_to_mrc("rel-"+ERRMAP_, err_rel, mrc)
# write output maps
print("\n Writing output density maps:")
print("%37s %s" % ("PLUMED output map :",PLUMED_))
print("%37s %s" % ("MRC output map :",PLUMED_.split(".")[0]+".mrc"))

# prepare empty data
data = np.zeros((mrc_p["nbin"][2], mrc_p["nbin"][1], mrc_p["nbin"][0]), dtype=np.float32)
for i in range(len(kcenters)):
    ii = kcenters[i]
    data[ii[2],ii[1],ii[0]] = dcenters[i]
# print MRC output file
write_data_to_mrc(PLUMED_.split(".")[0]+".mrc", data, mrc)

# print PLUMED output file
kfile=open(PLUMED_, "w")
# print header
kfile.write("%9s %6s  %13s %13s %13s      %10s     %10s\n" % ("#! FIELDS","Id","Pos_0","Pos_1","Pos_2","Density","Error"))
# cycle on voxels
for i in range(len(kcenters)):
    # voxel Id
    kfile.write("          %6d  " % i)
    # 3D index
    ii = kcenters[i]
    # voxel coordinates
    for j in range(3):
        # convert to nm
        xc = 0.1 * ( mrc_p["x0"][j] + float(ii[j]) * mrc_p["dx"][j] )
        kfile.write("%13.7lf " % xc)
    # density 
    ov = 1000.0 * dcenters[i]
    # print to file
    kfile.write("  %10.7e"  % ov)
    # error
    overr = 0.0
    if(do_halfmaps): overr = 1000.0 * err[ii[2],ii[1],ii[0]]
    kfile.write("  %10.7e\n" % overr)
print("\n All done!\n")
# close file
kfile.close()
