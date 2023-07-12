#!/bin/bash
#
# variables

#LOGFILE=/home/uttam/4D_recon/C/Patient_data/rest/2012041701_recon_rest_128x128x120_p0.log

# system matrix file=argv[1]
SMFILE=/home/uttam/4D_recon/C/hawkeye/even_no_pr_128x128x120-128x128x128.sm

#~/System_matrices/hawkeye/conct_no_pr_128x128x240-128x128x128.sm
TBFILE=./hawkeye/multihead_exponential_spline_GE.bin

# use the size of a single projection, i.e., 128x128 = 16384
SIZEOFASINGLEPROJECTIONANGLE=16384

NUMBEROFITERATIONS=54

# tomographic files =argv[2] #projection data from patient data
TOMOFILE=/home/uttam/4D_recon/C/Patient_data/2012041801/rest/2012041801_rest_128x128x120_p19.f32le

#f32le=float 32 little endian

# recon files=argv[3] #out put file
RECONFILE=/home/uttam/4D_recon/C/Patient_data/2012041801/rest/recon/even_2012041801_rest_128x128x120_p19.f32le

# recon

#echo ""
#echo "PPC dynamic recon of the 3600 frames of dynamic phantom data taken from 12 sets of even and odd rotations"
#echo ""

# for 3d put R 
# for 4d put r

#argv[4]=R (for 3d),r (for 4D)

#add these for 4D 
#argv[[6]
#$TBFILE
#argv[7] 
#$SIZEOFASINGLEPROJECTIONANGLE

echo "./use_sm $SMFILE $TOMOFILE $RECONFILE r $NUMBEROFITERATIONS $TBFILE $SIZEOFASINGLEPROJECTIONANGLE"

./use_sm $SMFILE $TOMOFILE $RECONFILE 	R $NUMBEROFITERATIONS 

#0	1	2	3   		4	5

echo ""

