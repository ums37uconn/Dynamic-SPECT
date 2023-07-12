#!/bin/bash
#

# system matrix file=argv[1]
#SMFILE=/home/uttam/4D_recon/C/ParallelPointResponseNew/system_matrix/2012031201_rest.sm

#SMFILE=/home/uttam/4D_recon/C/Phantom/mcat/system_matrix/phantom_amap.sm

SMFILE=/home/uttam/4D_recon/C/hawkeye/even_no_pr_128x128x120-128x128x128.sm

#TBFILE=/home/uttam/4D_recon/C/BSpline/phantom_basis_sine.bin

TBFILE=/home/uttam/4D_recon/C/BSpline/phantom_time_basis_1x15.bin

GBFILE=/home/uttam/4D_recon/C/BSpline/phantom_gate_basis_1x8.bin

# use the size of a single projection, i.e., 128x128 = 16384
SIZEOFASINGLEPROJECTIONANGLE=16384

NUMBEROFITERATIONS=1

# tomographic files =argv[2] #projection data from patient data
TOMOFILE=/home/uttam/4D_recon/C/Phantom/mcat/phantom_data/concatenated_8gate_1100-1290.bin

#f32le=float 32 little endian

# recon files=argv[3] #out put file
RECONFILE=/home/uttam/4D_recon/C/Phantom/mcat/Uttam/recon/concatenated_8gate_1100-1290.bin


# recon

#echo "./use_sm $SMFILE $TOMOFILE $RECONFILE r $NUMBEROFITERATIONS $TBFILE $SIZEOFASINGLEPROJECTIONANGLE"

#./use_sm $SMFILE $TOMOFILE $RECONFILE f $NUMBEROFITERATIONS $TBFILE $SIZEOFASINGLEPROJECTIONANGLE

./bin/use_sm_5D $SMFILE  $TOMOFILE  $RECONFILE g $NUMBEROFITERATIONS  $TBFILE $SIZEOFASINGLEPROJECTIONANGLE $GBFILE

#0	1	2	3   		4	5	6	7	9



