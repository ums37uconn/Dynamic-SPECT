#!/bin/bash
#

SMFILE=/data5/Yoonsuk/system_matrix/6th_simulation/5/phantom_120views.sm
# use the size of a single projection, i.e., 128x128 = 16384
SIZEOFASINGLEPROJECTIONANGLE=16384

NUMBEROFITERATIONS=50

TOMOFILE=/data5/Yoonsuk/6th_simulation/5_8to1_uniform_1.5cm/8to1_1.5cm_seq1-2_cycle1-8_Det0_DEW_128-128-2880_32le_flip_f.raw_120views

TBFILE=/data2/4D_recon/C/BSpline/xcat/singlehead_Bspline_9basis_Yoonsuk_rg_480.bin

#TBFILE=/data2/4D_recon/C/BSpline/time_1x6_respiratory_5x5_cardiac_12x12_sigma_4x8.bin

# recon files=argv[3] #out put file
RECONFILE=/data5/Yoonsuk/recon/6th_simulation/5/phantom_120views.f32le


./use_sm $SMFILE  $TOMOFILE  $RECONFILE r $NUMBEROFITERATIONS  $TBFILE $SIZEOFASINGLEPROJECTIONANGLE

#0	1	2	3   		4	5


#echo " "

