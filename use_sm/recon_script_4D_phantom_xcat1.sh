#!/bin/bash
#
for i in  0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39

do

SMFILE=/data5/Yoonsuk/system_matrix/6th_simulation/5/phantom$i.sm

# use the size of a single projection, i.e., 128x128 = 16384
SIZEOFASINGLEPROJECTIONANGLE=16384

NUMBEROFITERATIONS=100

TOMOFILE=/data5/Yoonsuk/6th_simulation/5_8to1_uniform_1.5cm/LM$i

#f32le=float 32 little endian

TBFILE=/data2/4D_recon/C/BSpline/xcat/singlehead_Bspline_9basis_Yoonsuk_rg_72.bin
#TBFILE=/data2/4D_recon/C/BSpline/time_1x6_respiratory_5x5_cardiac_12x12_sigma_4x8.bin

# recon files=argv[3] #out put file
RECONFILE=/data5/Yoonsuk/recon/6th_simulation/5/LM$i


./use_sm $SMFILE  $TOMOFILE  $RECONFILE r $NUMBEROFITERATIONS  $TBFILE $SIZEOFASINGLEPROJECTIONANGLE

#0	1	2	3   		4	5

done

#echo " "

