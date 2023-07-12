#!/bin/bash
#

ARG1=0 #StartTime
ARG2=480 #EndTime
ARG3=480 #NumberOfTimes


INPUT=/data5/Yoonsuk/recon/6th_simulation/5/phantom_120views.f32le.iter050

TBFILE=/data2/4D_recon/C/BSpline/xcat/singlehead_Bspline_9basis_Yoonsuk_rg_480.bin
	#/data2/4D_recon/C/BSpline/singlehead_Bspline_9basis_Yoonsuk_rg.bin

OUTPUT=/data5/Yoonsuk/recon/time_series/5/phantom_120views.f32le

./time_series_generator_xcat $ARG1 $ARG2  $ARG3  $INPUT $TBFILE $OUTPUT 


echo ""
