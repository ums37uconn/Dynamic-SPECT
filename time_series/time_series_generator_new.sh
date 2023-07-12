#!/bin/bash
#

ARG1=0 #StartTime
ARG2=1200 #EndTime
ARG3=1200 #NumberOfTimes

INPUT=/data5/Dynamic_Cardiac_SPECT/recon/stress_even_test8.f32le.iter054


#TBFILE=/home/uttam/4D_recon/C/BSpline/singlehead_my_spline_9basis_stress.bin

TBFILE=/data2/4D_recon/C/BSpline/singlehead_Bspline_9basis_1.bin

OUTPUT=/data5/Dynamic_Cardiac_SPECT/time_series/2019103101_stress_even_test8.f32le

#OUTPUT=/home/uttam/4D_recon/C/Patient_Data_Postprocessed/March_2012/2012032301_stress_128x128x120x10_even_newbasis.f32le

./time_series_generator_new $ARG1 $ARG2  $ARG3  $INPUT $TBFILE $OUTPUT 


echo ""
