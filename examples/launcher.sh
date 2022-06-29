#! /bin/bash

RUN_NAME=20220612.24x116_onlyP

MASKFILE=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/MASK/meshmask.nc

RSTBEFORE=/g100_scratch/userexternal/gbolzon0/ENSDAM/${RUN_NAME}/BEFORE/CHL_SUP/
RST_AFTER=/g100_scratch/userexternal/gbolzon0/ENSDAM/${RUN_NAME}/AFTER/CHL_SUP/
ENSDAM_INPUTDIR=/g100_scratch/userexternal/gbolzon0/${RUN_NAME}/ensdam_inputs/
OUTPUTDIR=/g100_scratch/userexternal/gbolzon0/ENSDAM/${RUN_NAME}/ensdam_outputs/

mkdir -p $ENSDAM_INPUTDIR $OUTPUTDIR
python netcdf_gen.py -m $MASKFILE -b $RSTBEFORE -a $RST_AFTER -o $ENSDAM_INPUTDIR -s 20190201 -e 20190331



for I in `ls $ENSDAM_INPUTDIR `; do
   ./example_from_netcdf.x $ENSDAM_INPUTDIR/$I > $OUTPUTDIR/$I
done
