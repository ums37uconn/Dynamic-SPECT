#!/bin/bash
#

ARG1=9 #kernelSize
ARG2=4 #FWHM
ARG3=95 #numberOfTimeFrames

INPUT=/data5/Dynamic_Cardiac_SPECT/time_series/time_series_2019071701_stress_128x128x128sx20_even25.f32le

OUTPUT=/data5/Dynamic_Cardiac_SPECT/time_series/filter/time_series_2019071701_stress_128x128x128sx20_even_4.f32le

./gaussianFilter4D $ARG1 $ARG2  $ARG3  $INPUT $OUTPUT 


#echo ""
