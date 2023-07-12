#!/bin/bash
#

# system matrix file=argv[1]
#SMFILE=/home/uttam/4D_recon/C/ParallelPointResponseNew/system_matrix/2012031201_rest.sm

SMFILE=/data5/Yoonsuk/system_matrix/test_rg3.sm

#~/System_matrices/hawkeye/conct_no_pr_128x128x240-128x128x128.sm
TBFILE=/data2/4D_recon/C/BSpline/singlehead_Bspline_9basis_Yoonsuk_rg.bin

#TBFILE=/home/uttam/4D_recon/C/BSpline/phantom_basis.bin

# use the size of a single projection, i.e., 128x128 = 16384
SIZEOFASINGLEPROJECTIONANGLE=16384

NUMBEROFITERATIONS=100

# tomographic files =argv[2] #projection data from patient data
#TOMOFILE=/data5/Yoonsuk/simulation3/seq_act2_cycle1-6_Det0_72GATE_40FRAME_128-128-2880_32le_flip_f.raw_rg0.f32le
TOMOFILE=/data5/Yoonsuk/simulation3/seq_act2_cycle1-6_Det0_72GATE_40FRAME_128-128-2880_32le_flip_f.raw_LM28

#f32le=float 32 little endian

# recon files=argv[3] #out put file
RECONFILE=/data5/Yoonsuk/recon/test_LM28.f32le


# recon

#./use_sm $SMFILE $TOMOFILE $RECONFILE f $NUMBEROFITERATIONS $TBFILE $SIZEOFASINGLEPROJECTIONANGLE
#3D
#./use_sm $SMFILE  $TOMOFILE  $RECONFILE R $NUMBEROFITERATIONS # $TBFILE $SIZEOFASINGLEPROJECTIONANGLE
#4D
./use_sm $SMFILE  $TOMOFILE  $RECONFILE r $NUMBEROFITERATIONS $TBFILE $SIZEOFASINGLEPROJECTIONANGLE

#0	1	2	3   		4	5

echo ""exit


