#!/bin/bash
#

SMFILE=/data5/Yoonsuk/system_matrix/test_beating_atn_1-40_128-128-128_32le_4.4_flip_rot.sm

# use the size of a single projection, i.e., 128x128 = 16384
SIZEOFASINGLEPROJECTIONANGLE=16384

NUMBEROFITERATIONS=54

TOMOFILE=/data5/Yoonsuk/5th_simulation_0.5deg_2.2mmPhantom/x3_seq1-4_cycle1-8_Det0_128-128-2880_32le_2.2_flip_f.raw

#f32le=float 32 little endian

#TBFILE=/data2/4D_recon/C/BSpline/xcat/time_1x6_respiratory_5x5_cardiac_12x12_sigma_4x8.bin
TBFILE=/data2/4D_recon/C/BSpline/time_1x6_respiratory_5x5_cardiac_12x12_sigma_4x8.bin

# recon files=argv[3] #out put file
RECONFILE=/data5/Yoonsuk/recon/test6D_new.raw


./use_sm $SMFILE  $TOMOFILE  $RECONFILE r $NUMBEROFITERATIONS  $TBFILE $SIZEOFASINGLEPROJECTIONANGLE

#0	1	2	3   		4	5

echo ""

