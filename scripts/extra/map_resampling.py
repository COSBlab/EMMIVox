import SimpleITK as sitk
import sys

##
## you need the following libraries:
##  conda install -c simpleitk simpleitk
##
## Written by Max Bonomi (mbonomi@pasteur.fr)
##

# load mrc
mrc = sitk.ReadImage(sys.argv[1], sitk.sitkFloat32)
# read resampling factor
res = int(sys.argv[2])

# extract info about spacing and number of bins
dx   = mrc.GetSpacing()
nbin = mrc.GetSize()

# define new spacing and number of bins
nbin_3D=[]; dx_3D=[]
for i in range(0,3):
    nbin_3D.append(res * nbin[i])
    dx_3D.append(dx[i] / float(res))

# create empty image for resampling
image_3D = sitk.Image(nbin_3D[0], nbin_3D[1], nbin_3D[2], sitk.sitkFloat32)
image_3D.SetOrigin(mrc.GetOrigin())
image_3D.SetSpacing(dx_3D)

# Resample and write
mrc_resampled = sitk.Resample(mrc, image_3D, sitk.Transform(3, sitk.sitkIdentity), sitk.sitkLinear, 0.0, mrc.GetPixelID())
sitk.WriteImage(mrc_resampled, 'resampled.mrc')
