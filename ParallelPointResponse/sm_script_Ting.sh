
#!/bin/bash
#
#
#input parameters file

INPUTFILE=/data2/4D_recon/C/ParallelPointResponseNew/parameter/test.txt

#input attentuation map

ATNMAP=/data5/Ting/20191023/4D_amap_rest/Ting_patient2_amap_rest_128x128x128.d32le

#output system matrix file
#OUTPUTFILE=/home/uttam/4D_recon/C/Phantom/mcat/system_matrix/phantom_amap.sm


OUTPUTFILE=/data5/Ting/20191023/system_matrix/test8.sm

./Parallel_ptresp  $OUTPUTFILE	$INPUTFILE    $ATNMAP

	#0		1	2		3   

#echo ""top
