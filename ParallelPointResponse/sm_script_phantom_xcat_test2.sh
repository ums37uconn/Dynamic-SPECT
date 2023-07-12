#!/bin/bash
#
#
#input parameters file

INPUTFILE=/data2/4D_recon/C/ParallelPointResponseNew/parameter/xcat/120_views.txt

#input attentuation map

ATNMAP=/data5/Yoonsuk/6th_simulation/5_8to1_uniform_1.5cm/map20.double


#output system matrix file
#OUTPUTFILE=/home/uttam/4D_recon/C/Phantom/mcat/system_matrix/phantom_amap.sm


OUTPUTFILE=/data5/Yoonsuk/system_matrix/6th_simulation/5/phantom_120views.sm


./Parallel_ptresp  $OUTPUTFILE	$INPUTFILE    $ATNMAP

#0		1	2		3


#echo ""top
