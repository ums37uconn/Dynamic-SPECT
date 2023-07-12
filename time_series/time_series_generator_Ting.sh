#!/bin/bash
#

ARG1=0 #StartTime
ARG2=600 #EndTime
ARG3=600 #NumberOfTimes

INPUT=/data5/Ting/20191023/recon/test1_singlehead_Ting_patient2_rest_128x128x120x20_even.f32le.iter050

#TBFILE=/home/uttam/4D_recon/C/BSpline/singlehead_my_spline_9basis_rest.bin

TBFILE=/data2/4D_recon/C/BSpline/singlehead_Bspline_9basis_Ting.bin

OUTPUT=/data5/Ting/20191023/recon/time_series/test1_singlehead_Ting_patient2_rest_128x128x120x20_even.f32le

#OUTPUT=/home/uttam/4D_recon/C/Patient_Data_Postprocessed/March_2012/2012032301_rest_128x128x120x10_even_newbasis.f32le

./time_series_generator $ARG1 $ARG2  $ARG3  $INPUT $TBFILE $OUTPUT 


echo ""
