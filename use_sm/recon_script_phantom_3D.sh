#!/bin/bash
#
# variables


SMFILE=/home/uttam/4D_recon/C/Phantom/MCat/system_matrix/360_views.sm

# use the size of a single projection, i.e., 128x128 = 16384
SIZEOFASINGLEPROJECTIONANGLE=16384

NUMBEROFITERATIONS=57

# tomographic files =argv[2] #projection data from patient data
TOMOFILE=/home/uttam/4D_recon/C/Phantom/MCat/Phantom_Data/respiratory/Projection/concatenate_360views.bin

#f32le=float 32 little endian

# recon files=argv[3] #out put file
RECONFILE=/home/uttam/4D_recon/C/Phantom/MCat/Phantom_Data/respiratory/Recon/3D_recon_concatenate_360views.f32le


# for 3d put R 
# for 4d put r

#argv[4]=R (for 3d),r (for 4D)

#add these for 4D 
#argv[[6]
#$TBFILE
#argv[7] 
#$SIZEOFASINGLEPROJECTIONANGLE

#echo "./use_sm $SMFILE $TOMOFILE $RECONFILE r $NUMBEROFITERATIONS $TBFILE $SIZEOFASINGLEPROJECTIONANGLE"

./use_sm $SMFILE $TOMOFILE $RECONFILE 	R  $NUMBEROFITERATIONS 

#0	1	2	3   		4	5

#echo ""

