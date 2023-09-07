# Single-structure and ensemble refinement with cryo-EM maps and EMMIVOX

This repository is organized in the following directories:
* `scripts`: python scripts used for preprocessing and analysis of EMMIVOX simulations
* `tutorials`: two complete tutorials for single-structure refinement and ensemble modelling
* `docs`: a few interesting (quite old) papers 
* `GROMACS`: old repository of GROMACS mdp files - not really needed

## **Software requirements**

 Make sure you have installed:

 * A modern C++ compiler that supports OpenMP parallelization and C++14. I am currently using gcc v.9.2.0. Tested also with gcc v.8.2.1;
 * MPI library/compilers for multi-replica ensemble simulations. I am currently using openmpi v.4.0.5. Tested also with v.3.1.6;
 * Cuda, needed by both GROMACS and PLUMED. I am currenly using v.11.1. The exact version depends a bit on how old your GPUs are. I am currently using NVIDIA A100;
 * LibTorch. I am currenly using v.1.8.0 downloaded [here](https://pytorch.org/get-started/locally/). Tested also with v.1.8.1.
   Make sure you download the C++ version (LibTorch, not pytorch) that is supported by the Cuda version you installed; 
 * XDR library, which you can install with [conda](https://anaconda.org/conda-forge/xdrfile) or compile youself from [source](ftp://ftp.gromacs.org/contrib/xdrfile-1.1.4.tar.gz). This is very useful in postprocessing as it allows PLUMED to write trajectories in GROMACS xtc format;
 * Conda to install the python libraries needed by the pre- and post-processing scripts. Have a look [here](https://gitlab.pasteur.fr/mbonomi/cryo-em-emmi/-/blob/master/scripts/README.md) for more info about the libraries that you need to install.
 * [Phenix](https://phenix-online.org/documentation/index.html) (any recent version), if you want to validate single-structure refinement. Not really needed for ensemble modelling.

## **PLUMED installation**

### 1. Getting PLUMED

Clone this private PLUMED version from Pasteur GitLab:

`git clone https://gitlab.pasteur.fr/mbonomi/plumed-indeep.git`

and switch to the isdb branch:

`git checkout isdb`

### 2. Configuring and hacking the `Makefile.conf`
 
First we configure PLUMED with:

`./configure --enable-xdrfile --enable-mpi`

Scroll the log and make sure that both MPI and XDR are found (you will see a warning if not).
Now comes the tricky part. We need to manually add to `Makefile.conf` the location of the libtorch library (plus a couple of minor edits). We are working to
make `configure` do this automatically, but we are not there yet. These are the modifications you need to do:
   * Substitute `-std=c++11` with `-std=c++14`;
   * Add `-D_GLIBCXX_USE_CXX14_ABI=1` to `CFLAGS`, `CXXFLAGS`, `CXXFLAGS_NOOPENMP`, and `LD`;
   * Add `-I/path_to_libtorch/include/torch/csrc/api/include/ -I/path_to_libtorch/include -I/path_to_libtorch/include/torch` to `CPPFLAGS` and `LD`;
   * Add `-L/path_to_libtorch/lib -ltorch -lc10 -lc10_cuda -ltorch_cpu -ltorch_cuda` to `DYNAMIC_LIBS`
   * Add `-Wl,-rpath,/path_to_libtorch/lib` to `DYNAMIC_LIBS` and `LD`

where `/path_to_libtorch/` is the root directory where libtorch is installed.

### 3. Compiling

Just type:

`make -j 16`

In this case, I am compiling in parallel using 16 CPU cores. You might want to change this number depending on your machine.
After succesfull compilation, you will need to source the `sourceme.sh` file to have PLUMED available in your shell, wherever you are. 

## **GROMACS installation**

### 1. Getting GROMACS

I am currenly using GROMACS v.2020.5, which you can download [here](https://manual.gromacs.org/2020.5/download.html).

### 2. Patching with PLUMED

Go to the root directory of GROMACS and patch with PLUMED:

`plumed patch -p --shared`

and select `gromacs-2020.5`.

### 3. Configuring and compiling

Follow GROMACS instructions to compile the code. This usually means creating a `build` directory in the GROMACS root directory and run `cmake` from within:

`cmake ../ -DCUDA_TOOLKIT_ROOT_DIR=/opt/gensoft/exe/cuda/11.1 -DGMX_GPU=ON -DCMAKE_INSTALL_PREFIX=/path/to/install/directory/gromacs-2020.5-plumed -DGMX_OPENMP=ON -DGMX_MPI=ON -DGMX_THREAD_MPI=OFF -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc`

Please note that you need to specify the directory where you have Cuda installed (`-DCUDA_TOOLKIT_ROOT_DIR`) and where you want to install GROMACS (`-DCMAKE_INSTALL_PREFIX`). 
Compile and install as usual:

`make -j 16`

`make install`

In this case, I am compiling in parallel using 16 CPU cores. You might want to change this number depending on your machine.

### 4. You are good to go

You will need to source this file in order to have GROMACS available in your shell, wherever you are:

`source /path/to/install/directory/gromacs-2020.5-plumed/bin/GMXRC`

## **Credits and contact**

This material has been prepared by Samuel Hoff and Max Bonomi (Institut Pasteur, France).
For any technical question, please write to mbonomi_at_pasteur.fr.
