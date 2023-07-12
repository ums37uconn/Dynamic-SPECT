#!/bin/bash
#
#input parameters file
#INPUTFILE=parms_phantom.txt
#INPUTFILE=parms_small_PP5.txt
#INPUTFILE=parms_small.txt

INPUTFILE=parm8.txt

#attenuation map
#ATNMAP=/home/uttam/ParallelPointResponseNew/cardiac_torso_amap_stress_64x64x80_fy_fz_rz90.f64le

#ATNMAP=/home/uttam/4D_recon/C/Patient_data/2012042401/stress/amap_stress/4D_amap_stress/2012042401_amap_stress_128x128x128.f32le

#ATNMAP=/home/uttam/4D_recon/C/hawkeye/amaps/2012052601_amap_stress_128x128x160_fy_fz_rz90.f64le

#ATNMAP=/home/uttam/4D_recon/C/Patient_data/2013022001/raw/amap_rotated/2013022001_amap_stress_128x128x128_fy_fz_rz90.f64le

#system matrix output

#OUTPUTFILE=/home/uttam/ParallelPointResponseNew/system_matrix/phantom-Max_Point49-atn_cutoff0.01-64x64x120-64x64x64.sm

OUTPUTFILE=/home/uttam/4D_recon/C/ParallelPointResponseNew/system_matrix/system_matrix_no_amap_parm8.sm

#echo "./use_sm $SMFILE $TOMOFILE $RECONFILE r $NUMBEROFITERATIONS $TBFILE $SIZEOFASINGLEPROJECTIONANGLE"

./Parallel_ptresp  $OUTPUTFILE	$INPUTFILE  #$ATNMAP

	#0		1		2	3   

#echo ""
