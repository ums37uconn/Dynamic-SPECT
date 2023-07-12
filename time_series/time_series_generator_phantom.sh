#!/bin/bash
#

ARG1=0 #StartTime
ARG2=1200 #EndTime
ARG3=1200 #NumberOfTimes

INPUT=/home/uttam/4D_recon/C/Phantom/MCat/Phantom_Data/respiratory/Recon/time_respiratory_singlehead_my_spline_20.bin.iter107
TBFILE=/home/uttam/4D_recon/C/BSpline/time_respiratory_singlehead_my_spline_20.bin
OUTPUT=/home/uttam/4D_recon/C/Phantom/MCat/Phantom_Data/respiratory/Recon/time_series/time_respiratory_singlehead_my_spline_20.f32le

./time_series_generator $ARG1 $ARG2  $ARG3  $INPUT $TBFILE $OUTPUT 


echo ""
