#!/bin/bash        
DDIR="../4-Production"
# extract NORM_DENSITY and RESOLUTION from `../1-Map-Preparation/log.preprocess` 
# and BEST_SCALE from ../3-Map-Scaling/BEST_SCALE
n=`grep NORM_DENSITY ../1-Map-Preparation/log.preprocess | awk '{print $NF}'`
r=`grep Resolution ../1-Map-Preparation/log.preprocess | awk '{print $NF/10.0}'`
s=`grep BEST_SCALE ../3-Map-Scaling/BEST_SCALE | awk '{print $NF}'`

# create plumed input file for production 
sed -e "s/NORM_DENSITY_/$n/g" plumed_EMMI_emin_template.dat | sed -e "s/RESOLUTION_/$r/g" | sed -e "s/SCALE_/$s/g" > plumed_EMMI_emin.dat

# Extract lowest energy frame from the single structure refinement (4-Production)
# get time (ps) of the frame with best score
b=`grep -v FIELDS ${DDIR}/COLVAR | sort -n -k 2 | head -n 1 | awk '{print $1}'`
# get the line of COLVAR corresponding to this frame
line=`awk '{print NR, $0}' ${DDIR}/COLVAR | grep -v FIELDS | sort -n -k 3 | head -n 1 | awk '{print $1}'`

# extract best frame (entire system)
echo 0 | gmx_mpi trjconv -f ${DDIR}/production.trr -o conf_best.gro -dump $b -s ${DDIR}/production.tpr

# create Bfactor file and take header
sed -n 1p ${DDIR}/EMMIStatus > EMMIStatus
# get the value of the line that has the lowest energy
sed -n ${line}p ${DDIR}/EMMIStatus >> EMMIStatus
