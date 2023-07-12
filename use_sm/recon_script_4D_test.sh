#!/bin/bash                                                                                                                                                    
#                                                                                                                                                              
#                                                                                                                                                              
#                                                                                                                                                              

SMFILE=/data5/Dynamic_Cardiac_SPECT/system_matrix/2019071701_rest.sm

TBFILE=/data2/4D_recon/C/BSpline/singlehead_my_spline_9basis.bin


SIZEOFASINGLEPROJECTIONANGLE=16384


NUMBEROFITERATIONS=1

TOMOFILE=/data5/Dynamic_Cardiac_SPECT/data/2019071701/4D_rest/2019071701_rest_128x128x120x10_even.f32le

RECONFILE=/home/uttam/scratch/output-recon.f32le


# for 3d put R                                                                                                                                                 
# for 4d put r                                                                                                                                                 

#argv[4]=R (for 3d),r (for 4D)                                                                                                                                 

#add these for 4D                                                                                                                                              
#argv[[6]                                                                                                                                                      
#$TBFILE                                                                                                                                                       
#argv[7]                                                                                                                                                       
#$SIZEOFASINGLEPROJECTIONANGLE                                                                                                                                 

#echo                                                                                                                                                          

./use_sm $SMFILE $TOMOFILE $RECONFILE r $NUMBEROFITERATIONS $TBFILE $SIZEOFASINGLEPROJECTIONANGLE


#./use_sm $SMFILE $TOMOFILE $RECONFILE R $NUMBEROFITERATIONS                                                                                                   


#0      1       2       3               4       5                                                                                                              

#echo "done" 
