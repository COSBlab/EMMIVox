#!/bin/bash
# extract NORM_DENSITY and RESOLUTION from `../1-Map-Preparation/log.preprocess` 
# and BEST_SCALE from ../3-Map-Scaling/BEST_SCALE
n=`grep NORM_DENSITY ../1-Map-Preparation/log.preprocess | awk '{print $NF}'`
r=`grep Resolution ../1-Map-Preparation/log.preprocess | awk '{print $NF/10.0}'`
s=`grep BEST_SCALE ../3-Map-Scaling/BEST_SCALE | awk '{print $NF}'`

# create plumed input file for production 
sed -e "s/NORM_DENSITY_/$n/g" plumed_EMMI_template.dat | sed -e "s/RESOLUTION_/$r/g" | sed -e "s/SCALE_/$s/g" > plumed_EMMI.dat
