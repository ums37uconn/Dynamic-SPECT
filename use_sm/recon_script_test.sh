#!/bin/bash
#
#
#


SMFILE=/data2/4D_recon/C/ParallelPointResponseNew/system_matrix/even/2012030501_amap_stress_128x128x128.sm
TOMOFILE=/data2/4D_recon/C/Patient_Data_Preprocessed/March_2012/0501/4D_stress/2012030501_stress_128x128x120x10_even.f32le
TBFILE=/data2/4D_recon/C/BSpline/singlehead_my_spline_9basis_stress.bin


SIZEOFASINGLEPROJECTIONANGLE=16384


NUMBEROFITERATIONS=1


RECONFILE=/data5/parametric/test


#strcpy(sysmat_fname, argv[1]);
#strcpy(sino_fname,   argv[2]);
#strcpy(image_fname,  argv[3]);
#strcpy(basis_fname,  argv[4]);
#Usage:./a.out sysmat.file sinogram.file image.file operation n_iterations {timebasis.file SingleProjectionSize} (Iteration numbers for output)
#Where "operation" can be either R for reconstruction, or F or B for forward or backprojection (CAPITAL LETTERS).

#echo

./use_sm $SMFILE $TOMOFILE $RECONFILE r  $NUMBEROFITERATIONS $TBFILE $SIZEOFASINGLEPROJECTIONANGLE


#./use_sm $SMFILE $TOMOFILE $RECONFILE 	R $NUMBEROFITERATIONS 


#0	1	2	3   		4	5

#echo ""

