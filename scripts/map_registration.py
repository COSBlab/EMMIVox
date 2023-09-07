import SimpleITK as sitk
import sys
import mrcfile
import numpy as np

##
## you need the following libraries:
##  conda install -c conda-forge mrcfile
##  conda install -c simpleitk simpleitk 
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

# load images
fixed_image =  sitk.ReadImage(sys.argv[1], sitk.sitkFloat32)
moving_map = sys.argv[2]
moving_image = sitk.ReadImage(moving_map, sitk.sitkFloat32)

# initialize
initial_transform = sitk.CenteredTransformInitializer(fixed_image, moving_image, sitk.Euler3DTransform(), sitk.CenteredTransformInitializerFilter.GEOMETRY)

# registration
registration_method = sitk.ImageRegistrationMethod()

# Similarity metric settings
registration_method.SetMetricAsCorrelation()
registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
registration_method.SetMetricSamplingPercentage(0.01)

# Interpolator
# they both perform well, but BSpline is much slower
interpolator = sitk.sitkGaussian
#interpolator = sitk.sitkBSpline
registration_method.SetInterpolator(interpolator)

# Optimizer settings
registration_method.SetOptimizerAsConjugateGradientLineSearch(learningRate=1.0, numberOfIterations=1000, convergenceMinimumValue=1e-6, convergenceWindowSize=10)

# Do the trick
registration_method.SetInitialTransform(initial_transform, inPlace=False)
final_transform = registration_method.Execute(fixed_image, moving_image)

# Always check the reason optimization terminated
print('Final metric value: {0}'.format(registration_method.GetMetricValue()))
print('Optimizer\'s stopping condition, {0}'.format(registration_method.GetOptimizerStopConditionDescription()))

# Resample and write
moving_resampled = sitk.Resample(moving_image, fixed_image, final_transform, interpolator, 0.0, moving_image.GetPixelID())
moving_map_al = '.'.join(moving_map.split(".")[:-1])+'-aligned.mrc'
sitk.WriteImage(moving_resampled, moving_map_al)

# now do some analysis with mrcfile
mrc_f = mrcfile.open(sys.argv[1], mode='r', permissive=True)
mrc_m = mrcfile.open('moved.mrc', mode='r', permissive=True)

# get parameters
mrc_f_p = get_map_parameters(mrc_f)
mrc_m_p = get_map_parameters(mrc_m)
# check parameters
for key in mrc_f_p:
    if(mrc_f_p[key]!=mrc_m_p[key]):
       print(" ERROR: map parameters are different!")
       exit()

# extract data and flatten array
data_f = mrc_f.data.flatten()
data_m = mrc_m.data.flatten()
# index of elements to keep
tokeep = np.where(data_f>float(sys.argv[3]))[0]
# delete the others
data_f = data_f[tokeep]
data_m = data_m[tokeep]
# define correlation (Chimera like)
cc = np.dot(data_f, data_m) / np.linalg.norm(data_f) / np.linalg.norm(data_m)
# printout
print("Correlation after registering: ",cc)
