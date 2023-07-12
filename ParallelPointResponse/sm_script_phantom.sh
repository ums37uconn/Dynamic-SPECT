#!/bin/bash
#
#
#input parameters file

INPUTFILE=/home/uttam/4D_recon/C/ParallelPointResponseNew/parameter/parms_forward.txt

#input attentuation map

ATNMAP=/home/uttam/4D_recon/C/Phantom/MCat/1061_atn_av_double.bin

#output system matrix file
#OUTPUTFILE=/home/uttam/4D_recon/C/Phantom/mcat/system_matrix/phantom_amap.sm



OUTPUTFILE=/home/uttam/4D_recon/C/Phantom/MCat/system_matrix/foward_projection/360.sm

./Parallel_ptresp  $OUTPUTFILE	$INPUTFILE    $ATNMAP

	#0		1	2		3   

#echo ""top
