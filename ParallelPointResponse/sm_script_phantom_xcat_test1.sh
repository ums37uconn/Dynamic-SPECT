#!/bin/bash
#
#
#input parameters file

for i in  0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	

do

INPUTFILE=/data2/4D_recon/C/ParallelPointResponseNew/parameter/xcat/angle_$i.txt

#input attentuation map

ATNMAP=/data5/Yoonsuk/6th_simulation/5_8to1_uniform_1.5cm/map$i.double


#output system matrix file
#OUTPUTFILE=/home/uttam/4D_recon/C/Phantom/mcat/system_matrix/phantom_amap.sm


OUTPUTFILE=/data5/Yoonsuk/system_matrix/6th_simulation/5/phantom$i.sm


./Parallel_ptresp  $OUTPUTFILE	$INPUTFILE    $ATNMAP

#0		1	2		3

done

#echo ""top
