# Single-structure model and ensemble refinement with cryo-EM maps and EMMIVox

Here you can find scripts and tutorials to perform single-structure model and ensemble refinement
using cryo-EM maps and the EMMIVox approach introduced in:

S. Hoff, F. E. Thomasen, K. Lindorff-Larsen, M. Bonomi. Accurate model and ensemble refinement using cryo-electron microscopy maps and Bayesian inference.
bioRxiv 2023 doi: [add ref](). 

This repository is organized in the following two directories:
* `scripts`: python scripts used for preprocessing and analysis of EMMIVox simulations
* `tutorials`: complete tutorials for single-structure and ensemble refinement

## **Software requirements**

 Make sure you have installed:

 * A modern C++ compiler that supports OpenMP parallelization and C++17.
 * MPI library/compilers for multi-replica ensemble simulations.
 * Cuda, needed by both GROMACS and PLUMED. The exact version depends a bit on how old your GPUs are.
 * [LibTorch](https://pytorch.org/get-started/locally/). Make sure you download the C++ version (LibTorch, not pytorch) that is supported by the Cuda version you installed; 
 * Conda to install the python libraries needed by the pre- and post-processing scripts. Have a look [here](https://github.com/maxbonomi/EMMIVox/tree/main/scripts) for more info about the libraries that you need to install.
 * [Phenix](https://phenix-online.org/documentation/index.html) (any recent version), if you want to validate single-structure refinement. Not really needed for ensemble modelling.

## **PLUMED installation**

### 1. Getting PLUMED

EMMIVox is currently implemented in the GitHub master branch of PLUMED, which you can obtain with:

`git clone https://github.com/plumed/plumed2.git`

or downloading the following zip archive:

`wget https://github.com/plumed/plumed2/archive/refs/heads/master.zip`

### 2. Configuring and compiling PLUMED
 
Please have a look [here](https://www.plumed.org/doc-master/user-doc/html/_i_s_d_b.html) for detailed instructions about compiling PLUMED with Libtorch support.
The main point is to enable Libtorch with:

`./configure --enable-libtorch`

## **GROMACS installation**

Detailed instructions about patching GROMACS with PLUMED, configuration and installation are available [here](https://www.plumed.org/doc-master/user-doc/html/_installation.html).

## **Credits and contact**

This material has been prepared by Samuel Hoff and Max Bonomi (Institut Pasteur, France).
For any technical question, please write to mbonomi_at_pasteur.fr.
