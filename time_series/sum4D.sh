#!/bin/bash
#
NUMBEROFTIMEFRAMES=95 #numberOfTimeFrames
#input 4D filtered data
FILTERDATA=/data5/Dynamic_Cardiac_SPECT/time_series/filter/2019071701_stress_128x128x128x95_Bspline_9basis_2_SC.f32le
#output 4D filtered data
OUTPUT=/data5/Dynamic_Cardiac_SPECT/time_series/filter/normalization/2019071701_stress_128x128x128x95_Bspline_9basis_2_SC.f32le
#output 3D sum data
OUTPUTSUM=/data5/Dynamic_Cardiac_SPECT/time_series/filter/sum4D/2019071701_stress_128x128x128x1_Bspline_9basis_2_SC.f32le

./sum4D $NUMBEROFTIMEFRAMES  $FILTERDATA $OUTPUT $OUTPUTSUM 


#echo ""
