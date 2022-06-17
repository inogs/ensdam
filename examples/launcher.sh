#! /bin/bash

INPUT_DIR=/g100_scratch/userexternal/gbolzon0/ENSDAM/20220612.24x116_onlyP/ensdam_inputs/
OUTPUTDIR=/g100_scratch/userexternal/gbolzon0/ENSDAM/20220612.24x116_onlyP/ensdam_outputs/

mkdir -p $OUTPUTDIR

for I in `ls $INPUT_DIR `; do
   ./example_from_netcdf.x $INPUT_DIR/$I > $OUTPUTDIR/$I
done
