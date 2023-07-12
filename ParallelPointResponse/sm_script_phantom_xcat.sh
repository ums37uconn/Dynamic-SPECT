#!/bin/bash
#
#
#input parameters file

INPUTFILE=/data2/4D_recon/C/ParallelPointResponseNew/parameter/parms_phantom.txt

#input attentuation map

ATNMAP=/data5/Yoonsuk/5th_simulation_0.5deg_2.2mmPhantom/beating_atn_1-40_128-128-128_32le_4.4_flip_rot.double

#output system matrix file
#OUTPUTFILE=/home/uttam/4D_recon/C/Phantom/mcat/system_matrix/phantom_amap.sm


OUTPUTFILE=/data5/Yoonsuk/system_matrix/test_beating_atn_1-40_128-128-128_32le_4.4_flip_rot.sm

./Parallel_ptresp  $OUTPUTFILE	$INPUTFILE    $ATNMAP

	#0		1	2		3   

#echo ""top
