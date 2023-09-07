# Python scripts for pre- and post-processing of EMMIVox simulations

Here you will find some python scripts useful to setup and analyze EMMIVox simulations. For info about their syntax, just type:

`python script.py -h`

Brief overview:
 * `align-VOXELS.py`: align two PDBs based on an atom selections, dump the transformation to file, apply transformation to a cryo-EM map, write transformed cryo-EM map to file;
 * `align-PDBs.py`: transform one PDB based on the inverse transformation read from file;
 * `align-XTC.py`: transform a trajectory in XTC format based on the inverse transformation read from file;
 * `map_registration.py`: align and resample a cryo-EM map on the grid of another cryo-EM map;
 * `cryo-EM_preprocess.py`: filter correlated voxels of a cryo-EM map, calculate error map from two half maps, zone the map close to the input model, write filtered map in PLUMED format to file;
 * `cryo-EM_validate.py`: calculate model/map fit (CCmask like) from one PDB, two PDBs, or two PDBs and a trajectory (average map);
 * `add-BFACT.py`: read a PDB file and a EMMIStatus file (with Bfactor) and create a new PDB file with Bfactor column.

## **Software installation**

You can create a conda environment to install all the python libraries needed to run pre- and post-processing tools.

* `conda create --name emmivox`

* `conda activate emmivox`

* `conda install -c conda-forge mrcfile mdanalysis biopython`

   Make sure you are installing MDAnalysis version >= 2.0.0.

  `conda install -c simpleitk simpleitk`

* You will also need pytorch for validating cryo-EM ensembles, not necessarily the same version as libtorch. You will find instructions [here](https://pytorch.org), the version probably depends on your Cuda installation.

  `conda install pytorch torchvision torchaudio cudatoolkit=11.3 -c pytorch`
