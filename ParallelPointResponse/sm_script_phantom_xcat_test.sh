#!/bin/bash
#
#
#input parameters file

for i in 341 342 343 344 345 346 347 348 349 350 351 352 353 354 355 356 357 358 359 360
      
do

INPUTFILE=/data2/4D_recon/C/ParallelPointResponseNew/parameter/xcat/angle_$i.txt

#input attentuation map

ATNMAP=/home/uttam/xcat/data/atn_crop_flip_64x64x64x16.d64le 

#output system matrix file
#OUTPUTFILE=/home/uttam/4D_recon/C/Phantom/mcat/system_matrix/phantom_amap.sm


OUTPUTFILE=/home/uttam/xcat/data/test_angle$i.sm

./Parallel_ptresp  $OUTPUTFILE	$INPUTFILE    $ATNMAP

#0		1	2		3

done

#echo ""top
