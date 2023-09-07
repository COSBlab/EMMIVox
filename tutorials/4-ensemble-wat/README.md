# Tutorial for ensemble modelling with EMMIVOX
These are the steps to do ensemble modelling using a cryo-EM map and EMMIVOX.
Each step of the procedure will be carried out in a separate directory.

**Note**: To do ensemble modelling you need to first complete the single-structure refinement tutorial.

## 0. Production 

   * Prepare the input files and directories for ensemble modelling with 16 replicas using:

     `bash prepare_PLUMED_Ensemble_input.sh 16`

     You can specify fewer or more replicas depending on the number of CPU cores and GPUs available.

     **Note**: If you have a monomeric protein or a heterocomplex, you need to edit `plumed_EMMI_template.dat` before executing the `prepare_PLUMED_Ensemble_input.sh` script
               and comment the line starting with `BFACT_NOCHAIN`. This option is used here since we are modelling 5 identical chains and
               we want the Bfactor of the same residue in different chains to be equal.

   * Run the simulation in parallel on the cluster, using 16 MPI processes, each one parallelized on multiple CPU cores (`$OMP_NUM_THREADS`).
     You might need to adapt this line depending on your command to submit parallel jobs, i.e. srun, mpiexec, or mpirun.

     `srun -N 16 gmx_mpi mdrun -pin on -ntomp $OMP_NUM_THREADS -plumed plumed_EMMI.dat -s production.tpr -multidir rep-*`

     **Note**: Usually a good setup for multiple-replica simulations is to allocate 1 GPU to each replica, and 6-10 CPU cores (`$OMP_NUM_THREADS`) per replica depending on the system size. You can safely allocate (using your job manager, i.e. slurm or pbs) 2 replicas on the same GPU without significant loss in performance.  

**Working directory**: `0-Production` 

## 1. Postprocessing and validation 

   * We first need to concatenate the trajectories of all replicas and fix discontinuities due to Periodic Boundary Conditions:

     `bash prepare_PLUMED_Ensemble_analysis.sh`

     Note: When prompted, choose the `System-MAP` group interactively to select the map atoms to be written in the output PDB file `conf_map.pdb`.

   * Now we need to transform back the concatenated trajectory `traj-all-PBC.xtc` to fit the original cryo-EM map 
     using the `transformation.dat` created during the map preparation stage of single-structure refinement:

     `python align-XTC.py conf_map.pdb traj-all-PBC.xtc traj-all-PBC-align.xtc ../../3-refinement-wat/1-Map-Preparation/transformation.dat` 

   * We add a Bfactor column to the reference PDB (`conf_map.pdb`). These Bfactors are the
     same for all residues, they were set before doing ensemble modelling equal to the minimum Bfactor found in the single-structure refinement:

     `python add-BFACT.py conf_map.pdb ../0-Production/EMMIStatus conf_map_bfact.pdb`

   * And finally we calculate model/map fit (Phenix CCmask-like) on the original PDB `7P6A.pdb` and on our ensemble:

     `python cryo-EM_validate.py ../../3-refinement-wat/1-Map-Preparation/emd_13223.map --pdbA=../../3-refinement-wat/1-Map-Preparation/7P6A.pdb --pdbC=conf_map_bfact.pdb --trjC=traj-all-PBC-align.xtc --threshold=0.0`

     **Note**: if you performed energy minimization and validation after single-structure refinement, you can also add this model to the comparison:

     `python cryo-EM_validate.py ../../3-refinement-wat/1-Map-Preparation/emd_13223.map --pdbA=../../3-refinement-wat/1-Map-Preparation/7P6A.pdb --pdbB=../../3-refinement-wat/5-Analysis/conf_map_bfact.pdb --pdbC=conf_map_bfact.pdb --trjC=traj-all-PBC-align.xtc --threshold=0.0`

**Working directory**: `1-Analysis`
