#!/bin/bash
#
#
#

#input parameter file
INPUTFILE=/data5/Dynamic_Cardiac_SPECT/parameter/parameter_2019071701_rest.txt

#attenuation map
#ATNMAP=/data2/4D_recon/C/Patient_Data_Preprocessed/Sept_Dec_2012/1210/4D_amap_rest/2012121001_amap_rest_128x128x128.d64le

ATNMAP=/data5/Dynamic_Cardiac_SPECT/data/2019071701/4D_amap_rest/2019071701_amap_rest_128x128x128.d64le

#system matrix output

OUTPUTFILE=/data5/Dynamic_Cardiac_SPECT/system_matrix/2019071701_rest_nohup.sm

#echo "./use_sm $SMFILE $TOMOFILE $RECONFILE r $NUMBEROFITERATIONS $TBFILE $SIZEOFASINGLEPROJECTIONANGLE"

./Parallel_ptresp  $OUTPUTFILE	$INPUTFILE  $ATNMAP

	#0		1		2	3   

#echo ""
