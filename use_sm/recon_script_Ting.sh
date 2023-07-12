#!/bin/bash
#

# system matrix file=argv[1]
#SMFILE=/home/uttam/4D_recon/C/ParallelPointResponseNew/system_matrix/2012031201_rest.sm

SMFILE=/data5/Ting/20191023/system_matrix/test5.sm

#~/System_matrices/hawkeye/conct_no_pr_128x128x240-128x128x128.sm
TBFILE=/data2/4D_recon/C/BSpline/singlehead_Bspline_9basis_Ting.bin

#TBFILE=/home/uttam/4D_recon/C/BSpline/phantom_basis.bin

# use the size of a single projection, i.e., 128x128 = 16384
SIZEOFASINGLEPROJECTIONANGLE=16384

NUMBEROFITERATIONS=50

# tomographic files =argv[2] #projection data from patient data
TOMOFILE=/data5/Ting/20191023/4D_rest/Ting_patient2_rest_128x128x120x20_odd.f32le

#f32le=float 32 little endian

# recon files=argv[3] #out put file
RECONFILE=/data5/Ting/20191023/recon/test5_singlehead_Ting_patient2_rest_128x128x120x20_odd.f32le


# recon

#./use_sm $SMFILE $TOMOFILE $RECONFILE f $NUMBEROFITERATIONS $TBFILE $SIZEOFASINGLEPROJECTIONANGLE
#3D
#./use_sm $SMFILE  $TOMOFILE  $RECONFILE R $NUMBEROFITERATIONS # $TBFILE $SIZEOFASINGLEPROJECTIONANGLE
#4D
./use_sm $SMFILE  $TOMOFILE  $RECONFILE r $NUMBEROFITERATIONS $TBFILE $SIZEOFASINGLEPROJECTIONANGLE

#0	1	2	3   		4	5

echo ""

