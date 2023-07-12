#!/bin/bash
#
#

SMFILE=/home/uttam/xcat/data/xcat_64x64x64x360.sm
RECONFILE=/home/uttam/data/concat_crop_flip_64x64x64x16.bin
TOMOFILE=/home/uttam/data/projection_test.bin
TBFILE=/data2/4D_recon/C/BSpline/phantom_basis001.bin
NUMBEROFITERATIONS=1

SINGLEPROJECTIONANGLE=1474560
#16384

#TBFILE=/home/uttam/4D_recon/C/BSpline/phantom_basis001.bin
# use the size of a single projection, i.e., 128x128 = 16384
#SIZEOFASINGLEPROJECTIONANGLE=16384

NUMBEROFITERATIONS=1

#tomographic files =argv[2] #projection data from patient data
#SMFILE=/home/uttam/4D_recon/C/Phantom/MCat/system_matrix/foward_projection/360.sm
#TOMOFILE=/home/uttam/4D_recon/C/Phantom/MCat/Phantom_Data/respiratory/10_act_60.bin
#RECONFILE=/home/uttam/4D_recon/C/Phantom/MCat/Phantom_Data/respiratory/Projection/10_act_60.bin
				   

# recon

#echo "./use_sm $SMFILE $TOMOFILE $RECONFILE r $NUMBEROFITERATIONS $TBFILE $SIZEOFASINGLEPROJECTIONANGLE"

#./use_sm $SMFILE $TOMOFILE $RECONFILE f $NUMBEROFITERATIONS $TBFILE $SIZEOFASINGLEPROJECTIONANGLE

./use_sm $SMFILE $TOMOFIL $RECONFILE f $NUMBEROFITERATIONS $TBFILE $SINGLEPROJECTIONANGLE

#0	1	2	3   		4	5

#echo ""

