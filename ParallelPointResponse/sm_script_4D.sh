#!/bin/bash
#
#
#

#input parameter file
INPUTFILE=/home/uttam/4D_recon/C/ParallelPointResponseNew/parameter/2012121001_stress.txt

#attenuation map
ATNMAP=/home/uttam/4D_recon/C/Patient_Data_Preprocessed/Sept_Dec_2012/1210/4D_amap_stress/2012121001_amap_stress_128x128x128.d64le

#system matrix output

OUTPUTFILE=/home/uttam/4D_recon/C/ParallelPointResponseNew/system_matrix/even/2012121001_amap_stress_128x128x128.sm

#echo "./use_sm $SMFILE $TOMOFILE $RECONFILE r $NUMBEROFITERATIONS $TBFILE $SIZEOFASINGLEPROJECTIONANGLE"

./Parallel_ptresp  $OUTPUTFILE	$INPUTFILE  $ATNMAP

	#0		1		2	3   

#echo ""
