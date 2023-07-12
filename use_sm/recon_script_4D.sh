#!/bin/bash
#
#
#


SMFILE=/home/uttam/4D_recon/C/ParallelPointResponseNew/system_matrix/even/2012032301_amap_stress_128x128x128.sm

TBFILE=/home/uttam/4D_recon/C/BSpline/singlehead_my_spline_9basis_stress.bin

SIZEOFASINGLEPROJECTIONANGLE=16384


NUMBEROFITERATIONS=54


TOMOFILE=/home/uttam/4D_recon/C/Patient_Data_Preprocessed/March_2012/2301/4D_stress/2012032301_stress_128x128x120x10_even.f32le

RECONFILE=/home/uttam/4D_recon/C/Patient_Data_Reconstructed/March_2012/2012032301_stress_128x128x120x10_even_newbasis.f32le


# for 3d put R 
# for 4d put r

#argv[4]=R (for 3d),r (for 4D)

#add these for 4D 
#argv[[6]
#$TBFILE
#argv[7] 
#$SIZEOFASINGLEPROJECTIONANGLE

#echo

./use_sm $SMFILE $TOMOFILE $RECONFILE r $NUMBEROFITERATIONS $TBFILE $SIZEOFASINGLEPROJECTIONANGLE


#./use_sm $SMFILE $TOMOFILE $RECONFILE 	R $NUMBEROFITERATIONS 


#0	1	2	3   		4	5

#echo ""

