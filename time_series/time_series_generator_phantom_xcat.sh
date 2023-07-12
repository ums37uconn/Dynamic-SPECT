#!/bin/bash
#

ARG1=0 #StartTime
ARG2=72 #EndTime
ARG3=72 #NumberOfTimes

for i in  0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39

do

INPUT=/data5/Yoonsuk/recon/6th_simulation/5/LM$((i)).iter100

TBFILE=/data2/4D_recon/C/BSpline/xcat/singlehead_Bspline_9basis_Yoonsuk_rg_72.bin
	#/data2/4D_recon/C/BSpline/singlehead_Bspline_9basis_Yoonsuk_rg.bin

OUTPUT=/data5/Yoonsuk/recon/time_series/5/LM$i

./time_series_generator_xcat $ARG1 $ARG2  $ARG3  $INPUT $TBFILE $OUTPUT 

echo $((i))

done

echo ""
