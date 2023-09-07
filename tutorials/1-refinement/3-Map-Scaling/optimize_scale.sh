#!/bin/bash
# number of CPU cores used by PLUMED 
export PLUMED_NUM_THREADS=$1
# extract NORM_DENSITY and RESOLUTION from `../1-Map-Preparation/log.preprocess` 
n=`grep NORM_DENSITY ../1-Map-Preparation/log.preprocess | awk '{print $NF}'`
r=`grep Resolution ../1-Map-Preparation/log.preprocess | awk '{print $NF/10.0}'`

# create a PDB file with only the XTC atoms
echo System-XTC | gmx_mpi trjconv -f ../2-Equilibration/em.gro -n ../0-Building/index.ndx -o step3_input_xtc.pdb -pbc nojump -s ../2-Equilibration/em.tpr

# loop over scale values
for d in $(seq 0.7 0.05 1.3)
do 
        # create directory and go into
	mkdir s-${d}; cd s-${d}
        # create plumed input file for postprocessing
        sed -e "s/NORM_DENSITY_/$n/g" ../plumed_EMMI_template_BFACT.dat | sed -e "s/RESOLUTION_/$r/g" | sed -e "s/SCALE_/$d/g" > plumed.dat
        # run PLUMED driver to calculate EMMIVOX score
        plumed driver --plumed plumed.dat --mf_xtc ../../2-Equilibration/nvt_posres.xtc > log.plumed
        # go back to root 	
	cd ../
done
# find scale value that minimize energy
cat s-*/COLVAR | grep -v FIELD | sort -n -k 2 | head -n 1 | awk '{print "BEST_SCALE: ",$3}' > BEST_SCALE
