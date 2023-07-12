#!/bin/bash
#

ARG1=0 #StartTime
ARG2=600 #EndTime
ARG3=600 #NumberOfTimes

INPUT=/home/uttam/scratch/Output4Drecon_patient2_even.f32le.iter010


#TBFILE=/home/uttam/4D_recon/C/BSpline/singlehead_my_spline_9basis_stress.bin

TBFILE=/data2/4D_recon/C/BSpline/singlehead_my_spline_9basis_stress.bin

OUTPUT=/home/uttam/scratch/time_series_Output4Drecon_patient2_even.f32le

#OUTPUT=/home/uttam/4D_recon/C/Patient_Data_Postprocessed/March_2012/2012032301_stress_128x128x120x10_even_newbasis.f32le

./time_series_generator $ARG1 $ARG2  $ARG3  $INPUT $TBFILE $OUTPUT 


echo ""
