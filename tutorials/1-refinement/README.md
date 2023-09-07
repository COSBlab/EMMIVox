# Tutorial for single-structure refinement with EMMIVox
These are the steps for refining a single structural model into a cryo-EM maps with EMMIVox.
Each step of the procedure will be carried out in a separate directory.

**Note**: all the python scripts are contained in [`scripts`](https://github.com/maxbonomi/EMMIVox/tree/main/scripts).

## 0. System setup

   * Create the topology files and the initial conformation in GROMACS format starting from the deposited PDB (`7P6A.pdb`). We have done this for you using [CHARMM-GUI](https://www.charmm-gui.org). CHARMM-GUI can also complete your system if there are missing residues and create the force field for a small-molecule, if present.

      **Note**: The conformation in the `step3_input.gro` file (GROMACS format) produced by CHARMM-GUI is identical to one in the deposited PDB.
                However, CHARMM-GUI will probably translate and rotate your initial PDB, which will then not fit the input cryo-EM map anymore. Furthermore, atoms order might be different.

   * CHARMM-GUI renumbers residues and chains in `step3_input.gro`, so we first need to make sure we have a PDB file consistent with the `gro`
     atom order. We can create this PDB file using this command:

     `bash renumber.sh step3_input.gro step3_input.pdb`

   * Add to the index file created by CHARMM-GUI (`index.ndx`) two custom groups:

     * `System-MAP`, which contains all the atoms that will be used to generate the cryo-EM map. You can do this with `make_ndx.py` in [`scripts`](https://github.com/maxbonomi/EMMIVox/tree/main/scripts) using [`MDAnalysis`](https://www.mdanalysis.org) selection syntax. Hydrogen atoms and the carboxylate oxygens of glutamic/aspartic acid will be automatically removed from this group, as they are not used in PLUMED to calculate the cryo-EM map. A second group, called `System-MAP-H` will also be created to include these missing atoms (mostly to write them in the trajectory file).

        `python make_ndx.py step3_input.gro "protein" System-MAP --ndx index.ndx`
       
        **Note**: if you have a small molecule or other non-protein (and non-water) atoms that you want to include in the EMMIVox restraint, please add them to the MDAnalysis selection.

     * `System-XTC`, which contains all the atoms that will be written to the GROMACS trajectory file (xtc format). In this case, we want to write all the atoms used for the cryo-EM calculation plus hydrogen and carboxylate oxygens of glutamic/aspartic acid.
 
        `python make_XTC_ndx.py System-MAP-H --ndx index.ndx`

   **Working directory**: `0-Building`

## 1. Map preparation

   * At this stage we need to download the cryo-EM full map [`emd_13223.map`](https://ftp.wwpdb.org/pub/emdb/structures/EMD-13223/map/emd_13223.map.gz),
     the two half-maps [`emd_13223_half_map_1.map`](https://ftp.wwpdb.org/pub/emdb/structures/EMD-13223/other/emd_13223_half_map_1.map.gz)
     and [`emd_13223_half_map_2.map`](https://ftp.wwpdb.org/pub/emdb/structures/EMD-13223/other/emd_13223_half_map_2.map.gz),
     and the PDB [`7P6A.pdb`](https://files.rcsb.org/download/7P6A.pdb), which will be needed to zone the density map close to the model (optional, but speeds up things a lot).

   * To convert the input cryo-EM map to PLUMED format, calculate an error map from the two half maps, and optionally filter voxels by correlation, you need to execute this:

      `python cryo-EM_preprocess.py emd_13223.map 0.9 emd_plumed.dat --zone 3.5 --zone_PDB 7P6A.pdb --zone_sel "protein" --halfmaps emd_13223_half_map_1.map emd_13223_half_map_2.map > log.preprocess` 

       **Note 1**: The two half maps should be registered (aligned and with the same grid parameters). If they are not, you can use this command to register them:

       `python map_registration.py emd.map emd_half_map_1.map`

       `python map_registration.py emd.map emd_half_map_2.map`

       **Note 2**: The value `0.9` is the cutoff to exclude correlated voxels (above this threshold). If you want to keep all the voxels of the input map, set this to `1.0`. `emd_plumed.dat` is the name of the output map in PLUMED format. `--zone` can be used to keep only voxels within a certain distance (here 3.5 Ang) from the model specified by `--zone_PDB`. 

       **Note 3**: If you want to give a look at the map after the preparation, you can inspect `emd_plumed.mrc`. This is the map that will be used in PLUMED for modelling and corresponds to `emd_plumed.dat`.

   * Now we align the cryo-EM map in PLUMED format (`emd_plumed.dat`) to the GROMACS conformation (`step3_input.pdb`). This is needed because CHARMM-GUI probably moved the initial structure during system setup. 

      `python align-VOXELS.py ../0-Building/step3_input.pdb 7P6A.pdb emd_plumed_aligned.dat emd_plumed.dat --ref_sel "backbone" --mobile_sel "backbone and not altLoc B and not (name O and resid 379)"`

       **Note 4**: Since there are some discrepancies between the number of atoms in the two PDBs, we have specified a selection of common atoms for alignment using [`MDAnalysis`](https://www.mdanalysis.org) selection syntax. `--ref_sel` refers to `step3_input.pdb` and `--mobile_sel` to `7P6A.pdb`. You can verify that the alignment is correct by checking the RMSD value in output: it should be zero as the two models are identical.

   **Working directory**: `1-Map-Preparation`

## 2. System equilibration

We need to prepare the system with an energy minimization and equilibration at room temperature. No cryo-EM restraints will be used at this stage.

   * Run energy minimization:

     `gmx_mpi grompp -f 0-em-steep.mdp -c ../0-Building/step3_input.gro -p ../0-Building/topol.top -o em.tpr`

     `gmx_mpi mdrun -pin on -deffnm em`

   * Do a 1-ns NPT equilibration with positional restraints on everything except water and ions:

     `gmx_mpi grompp -f 1_npt_posres.mdp -c em.gro  -n ../0-Building/index.ndx -p ../0-Building/topol.top -r em.gro -o npt_posres.tpr -maxwarn 1`

     `gmx_mpi mdrun -pin on -deffnm npt_posres -nsteps 500000`

   * Do a 2-ns NVT equilibration with positional restraints on everything except water and ions:

     `gmx_mpi grompp -f 2_nvt_posres.mdp -c npt_posres.gro -n ../0-Building/index.ndx -p ../0-Building/topol.top -r em.gro -o nvt_posres.tpr -maxwarn 1`

     `gmx_mpi mdrun -pin on -deffnm nvt_posres -nsteps 1000000`

   **Working directory**: `2-Equilibration`

## 3. Optimizing the cryo-EM map scaling factor

   We postprocess the NVT equilibration trajectory to obtain an optimal scaling factor between experimental and predicted map. 
   
   `bash optimize_scale.sh 8`

   Here we are using 8 CPU cores (and the GPU) to postprocess the trajectory with PLUMED. The optimal value of the scale will be printed in `BEST_SCALE` at the end of the optimization. This will be used for both single-structure refinement and ensemble modelling.

   **Note**: If you have a monomeric protein or a heterocomplex, you need to edit `plumed_EMMI_template_BFACT.dat` before executing the `optimize_scale.sh` script
               and comment the line starting with `BFACT_NOCHAIN`. This option is used here since we are modelling 5 identical chains and 
               we want the Bfactor of the same residue in different chains to be equal.

   **Working directory**: `3-Map-Scaling`

## 4. Production

   To run the production simulation for single-structure refinement we need to:

   * Prepare the `plumed_EMMI.dat` file with production parameters using:

     `bash prepare_PLUMED_input.sh`

     **Note**: If you have a monomeric protein or a heterocomplex, you need to edit `plumed_EMMI_template.dat` before executing the `prepare_PLUMED_input.sh` script
               and comment the line starting with `BFACT_NOCHAIN`. This option is used here since we are modelling 5 identical chains and
               we want the Bfactor of the same residue in different chains to be equal.

   * Run a 10-ns long production run following the instructions below, after setting the number of CPU cores to use (`$OMP_NUM_THREADS`). 

     `gmx_mpi grompp -f 4-nvt-production.mdp -c ../2-Equilibration/nvt_posres.gro -n ../0-Building/index.ndx -p ../0-Building/topol.top -o production.tpr`

     `gmx_mpi mdrun -pin on -deffnm production -ntomp $OMP_NUM_THREADS -nsteps 5000000 -plumed plumed_EMMI.dat`

   **Working directory**: `4-Production`

## 5. Postprocessing and validation

   **Note**: You do not need to perfom these actions if you are interested only in ensemble modelling.
    
   * We first need to extract the conformation with best EMMIVox score from our production run and perform a short energy minimization. To setup all the files needed for minimization, you need to execute this script:
     
     `bash prepare_PLUMED_input_emin.sh`

     **Note**: If you have a monomeric protein or a heterocomplex, you need to edit `plumed_EMMI_emin_template.dat` before executing the `prepare_PLUMED_input_emin.sh` script
               and comment the line starting with `BFACT_NOCHAIN`. This option is used here since we are modelling 5 identical chains and
               we want the Bfactor of the same residue in different chains to be equal. 

   * Now we run energy minimization using 8 CPU cores:
    
     `bash run_PLUMED_emin.sh 8`

   * After minimization is complete, we need to align `conf_pbc.pdb` to the original cryo-EM map using the `transformation.dat`
     file created during the map preparation stage:

     `python align-PDBs.py conf_pbc.pdb conf_pbc_aligned.pdb ../1-Map-Preparation/transformation.dat`

   * Finally, we add the BFactors column to our aligned model:

     `python add-BFACT.py conf_pbc_aligned.pdb EMMIStatus conf_phenix.pdb`

   * The output PDB file `conf_phenix.pdb` is ready to be validated with PHENIX:

     `bash do_PHENIX conf_phenix.pdb ../1-Map-Preparation/emd_13223.map 1.9 > results.PLUMED`
   
     where `1.9` is the resolution of the input map `emd_13223.map` in Angstrom. Validation metrics are saved in `results.PLUMED`.

   * Compare the EMMIVOX-refined model with the deposited PDB:

      `bash do_PHENIX ../1-Map-Preparation/7P6A.pdb ../1-Map-Preparation/emd_13223.map 1.9 > results.PDB`

   **Working directory**: `5-Analysis`
