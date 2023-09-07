# 1) concatenate all the trajectories from multiple replicas -> traj-all.xtc
gmx_mpi trjcat -f ../0-Production/rep-*/traj*.xtc -cat -o traj-all.xtc

# 2) fix PBCs with PLUMED and select only map atoms -> traj-all-PBC.xtc
plumed driver --plumed plumed_fix_pbc.dat --mf_xtc traj-all.xtc

# 3) extract from step3_input_xtc.pdb only the atoms used to generate the cryo-em map. 
gmx_mpi trjconv -f ../../3-refinement-wat/3-Map-Scaling/step3_input_xtc.pdb -s ../0-Production/rep-00/production.tpr -n ../../3-refinement-wat/0-Building/index.ndx -o conf_map.pdb
