#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


// input sinogram can be either 16-bit unsigned integer (unsigned short) or float
#define FLOAT_INPUT_SINOGRAM
// uncomment this, if input sinogram is float, leave like this otherwise

// Volume is best output as 4-bit float
#define FLOAT_INPUT_VOLUME

// Max iterations and outputs of different iterations
#define MAX_ITERATIONS 1000
#define MAX_OUTPUTS 100

// Max number of rotations of camera accepted by the forward/backward projection with a time basis
#define	MAX_ROTATIONS	1000

// This defines how small zero is, and what is one divided by zero
#define ONE_OVER_EPSILON 1.0e+13
#define EPSILON 1.0e-13

//.......for dynamic recon (4D)

	void TB_MLEM_recon(int Nsinogram, int Nvolume, int singleProjectionSize, 
			int basisNum, int basisLen, float *sinogram, float *coeffVolume, 
			float **basis, char *sysmat_fname, int LastIteration, int OutputNumber, 
			int *IterationsForOutput, char *image_fname);
	void TimeBasisBP(int Nsinogram, int Nvolume, int singleProjectionSize, 
			int basisNum, int basisLen, float *sinogram, float *coeffVolume, 
			float **basis, char *sysmat_fname);
	void TimeBasisFP(int Nsinogram, int Nvolume, int singleProjectionSize, int basisNum, 
			int basisLen, float *sinogram, float *coeffVolume, float **basis, 
			char *sysmat_fname);

//........for gate + dynamic recon (5D)


	void GB_MLEM_recon(int Nsinogram, int Nvolume, int singleProjectionSize, 
			int basisNum, int basisLen,int gatebasisNum, int gatebasisLen, float *sinogram, float *coeffVolume, 
			float **basis,float **gatebasis, char *sysmat_fname, int LastIteration, int OutputNumber, int *IterationsForOutput, char *image_fname);

	void GateBasisBP(int Nsinogram, int Nvolume, int singleProjectionSize, int basisNum, 
			int basisLen, int gatebasisNum, int gatebasisLen, float *sinogram, float *coeffVolume, float **basis, float **gatebasis,
			char *sysmat_fname);

	void GateBasisFP(int Nsinogram, int Nvolume, int singleProjectionSize, int basisNum, 
			int basisLen, int gatebasisNum, int gatebasisLen, float *sinogram, float *coeffVolume, float **basis, float **gatebasis,
			char *sysmat_fname);


int main(int argc, char **argv){
  
		int n, Nvolume, Nsinogram;
  		int IterationsForOutput[MAX_OUTPUTS];
  		int OutputNumber=1;
  		int index;

 		int tbFlag, tbNum, tbLen, singleProjectionSize;
		int gbFlag, gbNum, gbLen;
  		
		FILE *fid;
  		float **basis, **gate_basis;
  		float *sinogram, *volume;

  		char sino_fname[128], image_fname[128], sysmat_fname[128], basis_fname[128], gate_basis_fname[128], Action;
  		unsigned short tmp;

  		int minNumberOfArguments = 6;
 		int maxNumberOfArguments = 9;

  		if((argc < minNumberOfArguments ) || (argc > maxNumberOfArguments + MAX_OUTPUTS)){
   			 fprintf(stderr, "\nUsage:\n%s sysmat.file sinogram.file image.file" 
				"operation n_iterations {timebasis.file SingleProjectionSize}"
 				"(Iteration numbers for output)\n",argv[0]);
			
   		 	 exit(1);
 		 	}

 		strcpy(sysmat_fname, argv[1]);
  		strcpy(sino_fname,   argv[2]);
 		strcpy(image_fname,  argv[3]);
 		
		Action = argv[4][0];
		if(Action == 'r'){tbFlag = 1; gbFlag = 0;}
			else if(Action == 'g'){tbFlag = 1; gbFlag = 1;}
				else {
					fprintf(stderr, "unrecongnized option \"%s\"\n", argv[4]);
					exit(1);
				     }

 	 // get number of iterations to compute
 		int LastIteration = atoi( argv[minNumberOfArguments-1] );

	//open system matrix file
  	// fprintf( stderr, "Reading sinogram and image size from system matrix file...\n" );
 	 	if (( fid = fopen(sysmat_fname, "rb" )) == NULL ) {
      			fprintf(stderr, "Could not open sysmat file \"%s\" for reading\n", sysmat_fname);
      			exit(1);
  		}

  		fread( &Nsinogram, sizeof( int ), 1, fid );
  		fread( &Nvolume, sizeof( int ), 1, fid );
  		fclose( fid );
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Read time basis file
  	if ( tbFlag==1 ){ 
		  	strcpy( basis_fname,  argv[ minNumberOfArguments ] );
		 	if((fid = fopen(basis_fname, "rb"))==NULL)
				{
			 	 fprintf(stderr, "Could not open basis file %s for reading.\n",basis_fname);
			 	 exit(1);
		  		}

		  	fread(&tbNum, sizeof(int), 1, fid);
		 	fread(&tbLen, sizeof(int), 1, fid);

		 	 basis = (float **)   malloc(sizeof(float *)  *  tbNum);

		 	 for(index=0;index<tbNum;index++)
			 	 basis[index] = (float *) malloc(sizeof(float)  *  tbLen);

		 	 for(index=0;index<tbNum;index++)
			  	if(fread(&basis[index][0], sizeof(float), tbLen, fid)!=tbLen)
		 	  		{
					fprintf(stderr, "Basis file %s is corrupt, not long enough.\n",basis_fname);
					exit(1);
			  		}
		  	fclose(fid);
		 	singleProjectionSize = atoi( argv[ minNumberOfArguments + 1 ] );

		  	if((singleProjectionSize<=0)||(singleProjectionSize>Nsinogram))
				{
				fprintf(stderr, "SingleProjectionSize %d must be an integer smaller than the" 
					"size of the original sinogram %d.\n",singleProjectionSize,Nsinogram);
					exit(1);
		  		}

	  	Nsinogram = singleProjectionSize * tbLen;
	  	Nvolume = Nvolume * tbNum;

	//fprintf(stderr, "Nsinogram = %d\n",Nsinogram);

	}//tbFlag loop closes here

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//If 5D then read also gate basis file
  	if ( gbFlag==1 ){ 
			
			strcpy( gate_basis_fname,  argv[ 8 ] );//argument 8 is gate basis

		 	if((fid = fopen(gate_basis_fname, "rb"))==NULL)
				{
			 	 fprintf(stderr, "Could not open gate basis file %s for reading.\n",gate_basis_fname);
			 	 exit(1);
		  		}

			fread(&gbNum, sizeof(int), 1, fid);
		 	fread(&gbLen, sizeof(int), 1, fid);

		 	gate_basis = (float **)  malloc(sizeof(float *)  *  gbNum);

		 	for(index=0;index<gbNum;index++)
			 	 gate_basis[index] = (float *) malloc(sizeof(float)  *  gbLen);

		 	for(index=0;index<gbNum;index++)
			  	if(fread(&basis[index][0], sizeof(float), gbLen, fid)!=gbLen)
		 	  		{
						fprintf(stderr, "gate Basis file %s is corrupt, not long enough.\n",gate_basis_fname);
						exit(1);
			  		}
		  	fclose(fid);

		 	singleProjectionSize = atoi( argv[ minNumberOfArguments + 1 ] );

		  	if((singleProjectionSize<=0)||(singleProjectionSize>Nsinogram))
				{
				fprintf(stderr, "SingleProjectionSize %d must be an integer smaller than the" 
					"size of the original sinogram %d.\n",singleProjectionSize,Nsinogram);
					exit(1);
		  		}

	  	Nsinogram = singleProjectionSize * tbLen * gbLen;
	  	Nvolume = Nvolume * tbNum * gbNum;

	fprintf(stderr, "tbLen = %d  gbLen = %d \n", tbLen, gbLen);
	fprintf(stderr, "tbNum = %d  gbNum = %d \n", tbNum, gbNum);
	fprintf(stderr, "Nsinogram = %d\n",Nsinogram);

	} //gbFlag loop closes here


////////////////////////////////////////////////////////////////////////////////////////////

  		
	  		IterationsForOutput[0]=LastIteration;
		    

///////////////////////////////////////////////////////////////////////////////////////////////////////////
  
		sinogram = (float *) malloc(sizeof(float)*Nsinogram);
	  	volume   = (float *) malloc(sizeof(float)*Nvolume);


 		//read sinogram
    		if((fid = fopen(sino_fname, "rb"))==NULL){ 
			fprintf(stderr, "Could not open sinogram data file \"%s\" for reading\n", sino_fname);
     			exit(1);
			}
			#ifndef FLOAT_INPUT_SINOGRAM
    			for(n=0; n<Nsinogram; n++){ 			
				fread(&tmp, 2, 1, fid); // 2 = sizeof(unsigned short) == 16 bit
     				sinogram[n] = tmp; 
				}
			#else
   			fread(sinogram, sizeof(float), Nsinogram, fid);
			#endif
   		 fclose(fid);

////////////////////////////////////////////////////////////////////////////////////////////////////

// At this point, we have read all the required input data, except sysmat

  	
			if(Action == 'r')
	      			TB_MLEM_recon(Nsinogram, Nvolume, singleProjectionSize, 
					tbNum, tbLen, sinogram, volume, basis, sysmat_fname, 
					LastIteration, OutputNumber, IterationsForOutput, image_fname);


			if(Action == 'g')
	      			GB_MLEM_recon(Nsinogram, Nvolume, singleProjectionSize, 
					tbNum, tbLen,gbNum, gbLen, sinogram, volume, basis,gate_basis, sysmat_fname, 
					LastIteration, OutputNumber, IterationsForOutput, image_fname);
			

  	//free memory and return
  			free(volume);
  			free(sinogram);

  
  	if(tbFlag==1){
	 	for(index=0;index<tbNum;index++)
			 free(basis[index]);
			 free(basis);
 		 }

	if(gbFlag==1){
	 	for(index=0;index<gbNum;index++)
			 free(gate_basis[index]);
			 free(gate_basis);
 		 }

 	 return 1;

}
/////end of main/////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////
void GB_MLEM_recon(int Nsinogram, int Nvolume, int singleProjectionSize, 
			int basisNum, int basisLen,int gatebasisNum, int gatebasisLen, float *sinogram, float *coeffVolume, 
			float **basis,float **gatebasis, char *sysmat_fname, int LastIteration, int OutputNumber, 
			int *IterationsForOutput, char *image_fname){
 

		FILE *fid;
  		int n, niter,index;
  		float *tmp_volume, *tmp_sinogram, *constant_denominator;
  		char temp_fname[128];
  		// --- initial volume assignment: all pixels are one
  		for(n=0; n<Nvolume; n++) coeffVolume[n] = 1.0;

		// ------ create ml_em variables:
  		tmp_sinogram = (float *) malloc((size_t) (sizeof(float) * Nsinogram));
  		constant_denominator = (float *) malloc((size_t) (sizeof(float) * Nvolume));
  		tmp_volume = (float *) malloc((size_t) (sizeof(float) * Nvolume));

		// --- compute  element-by-element inverse of efficiency matrix
  		for(n=0; n<Nsinogram; n++) tmp_sinogram[n] = 1.0;

  		GateBasisBP(Nsinogram, Nvolume, singleProjectionSize, basisNum,
				basisLen, gatebasisNum, gatebasisLen, tmp_sinogram, coeffVolume, 
				basis, gatebasis, sysmat_fname);

/*
  		for(n=0; n<Nvolume; n++)
    		if(constant_denominator[n] > EPSILON)  constant_denominator[n] = 1./ constant_denominator[n];
    		else constant_denominator[n] = ONE_OVER_EPSILON;


	//  -------- ITERATION LOOP --------
  		for(niter=1; niter<=LastIteration; niter++){

    			//fprintf(stderr, "Iteration No %d of %d\n", niter, LastIteration);
			// compute the reprojection through the n-1 version of the file into tmp_sinogram

    		GateBasisFP(Nsinogram, Nvolume, singleProjectionSize, basisNum,
		basisLen, gatebasisNum, gatebasisLen, tmp_sinogram, coeffVolume, basis, gatebasis, sysmat_fname);

	// divide the sinogram by the tmp_sinogram
    	for(n=0; n<Nsinogram; n++){
      		if(sinogram[n] <= 0.0)  tmp_sinogram[n] = 0.0;
        	else if(tmp_sinogram[n] > EPSILON) tmp_sinogram[n] = sinogram[n] / tmp_sinogram[n];
           	else                          tmp_sinogram[n] = sinogram[n] * ONE_OVER_EPSILON;
		}


		// backproject the result into tmp_volume
      		GateBasisBP(Nsinogram, Nvolume, singleProjectionSize, basisNum,
				 basisLen, gatebasisNum,
				 gatebasisLen, tmp_sinogram, coeffVolume, basis, gatebasis, sysmat_fname);

		// multiply by the constant denominator
   		 for(n=0; n<Nvolume; n++)
    			{
      			coeffVolume[n] *=  constant_denominator[n] * tmp_volume[n];
      			if(coeffVolume[n] < 0.) coeffVolume[n] = 0.;
    			}


    		for(index=0;index<OutputNumber;index++)
			if(niter==IterationsForOutput[index])
    			{
				sprintf(temp_fname,"%s.iter%03d",image_fname,niter);
				if((fid = fopen(temp_fname, "wb"))==NULL)  //sprintf("%s.iter%d",image_fname,niter)
					{fprintf(stderr, "Could not open image file \"%s\" for writing\n", temp_fname);
					exit(1);}
			fwrite(coeffVolume, sizeof(float), Nvolume, fid);
			fclose(fid);
			}


  	}
	//end iteration here

*/
  	// end: free up memory
 	 free(tmp_sinogram);
 	 free(tmp_volume);
  	 free(constant_denominator);

 
		
}


//////////////// Backward projection using a system-matrix, modified by a set of time basis functions /////////////////

void GateBasisBP(int Nsinogram, int Nvolume, int singleProjectionSize, int basisNum, 
			int basisLen, int gatebasisNum, int gatebasisLen, float *sinogram, float *coeffVolume, float **basis, float **gatebasis,
			char *sysmat_fname){


  		int orgNvolume;
  		int ns, n, Nchunk, basisIndex, rotationIndex, rotationNum;
  		int currentProjection, placeInCurrentProjection, orgNsinogram, singleProjInRotation;
  		float *SMtemp;
  		int *Itemp,*basisInVolume,*rotationInSinogram;
  		int **singleProjFlag, **singleProjInBasis;
  		FILE *fid;
  		size_t size_int,size_intptr, size_float;

  		size_int = sizeof(int);
  		size_intptr = sizeof(int *);
  		size_float = sizeof(float);

  		// Check if paramters match for the total sinogram size 
  		if(singleProjectionSize*basisLen*gatebasisLen!=Nsinogram)
  			{
    				fprintf(stderr, "Number of single projections %d and time frames %d does not match output sinogram size %d\n", 							singleProjectionSize,basisLen,Nsinogram);
    				exit(-1);
  			}


		//open system matrix file
  		if((fid = fopen(sysmat_fname, "rb"))==NULL)
  			{
    				fprintf(stderr, "Could not open sysmat file %s\n", sysmat_fname);
    				exit(-1);
  			}

  	//Read in static sinogram size from original system matrix 
  		fread(&orgNsinogram, size_int, 1, fid);

  		singleProjInRotation=0;
  		while(singleProjInRotation*singleProjectionSize<orgNsinogram)
	  		singleProjInRotation++;


		fprintf(stderr, "orgNsinogram = %d\n",orgNsinogram);

  			currentProjection=0;
 			rotationNum=0;
  			while(currentProjection<basisLen*gatebasisLen){
	  				currentProjection+=singleProjInRotation;
	  				rotationNum++;
  				}
			fprintf(stderr, "currentProjection = %d\n",currentProjection);
			fprintf(stderr, "rotationNum = %d\n",rotationNum);


   		if(rotationNum>=MAX_ROTATIONS){
    			fprintf(stderr, "Using original sinogram size %d (from system matrix) for basis length %d exceeds %d rotatioms.\n", 
			orgNsinogram, basisLen, MAX_ROTATIONS);
    			exit(-1);
  			}
 		
  		fread(&orgNvolume, size_int, 1, fid);
		fprintf(stderr, "orgNvolume = %d\n",orgNvolume);

  		if(Nvolume!=orgNvolume*basisNum*gatebasisNum)
 		 {
    			fprintf(stderr, "Given Nvolume %d is not equal to expected %d * %d = %d.\n", Nvolume, orgNvolume, basisNum, orgNvolume*basisNum);
    			exit(-1);
  		 }
		
  		memset(coeffVolume, 0, size_float*Nvolume);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  		// Allocate memory and compute variables to speed up calculations
  		Itemp = (int *)   malloc(size_int  *  orgNvolume);
 		SMtemp = (float *) malloc(size_float * orgNvolume);

  		basisInVolume = (int *)   malloc(size_int  *  basisNum);

  		for(basisIndex=0;basisIndex<basisNum;basisIndex++)
	  		basisInVolume[basisIndex]=basisIndex*orgNvolume;

  		rotationInSinogram = (int *)   malloc(size_int  *  rotationNum);

  		for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++)
	  		rotationInSinogram[rotationIndex]=rotationIndex*orgNsinogram;

  			singleProjInBasis = (int **)   malloc(size_intptr  *  rotationNum);
  			singleProjFlag = (int **)   malloc(size_intptr  *  rotationNum);


  		for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++){
	  		singleProjInBasis[rotationIndex] = (int *)   malloc(size_int  *  singleProjInRotation);
	  		singleProjFlag[rotationIndex] = (int *)   malloc(size_int  *  singleProjInRotation);

	  		for(n=0;n<singleProjInRotation;n++){	  		
				singleProjInBasis[rotationIndex][n] = rotationIndex*singleProjInRotation + n;
				if(singleProjInBasis[rotationIndex][n]<basisLen)
					singleProjFlag[rotationIndex][n] = 1;
				else
				singleProjFlag[rotationIndex][n] = 0;
	  			}
 			 }

  			currentProjection=0;
  			placeInCurrentProjection=0;
 		
  		for(ns=0;ns<orgNsinogram;ns++){  		
   		 	fread(&n, size_int, 1, fid);
   		 		if(n!=ns){
      					fprintf(stderr, "Read in sinogram index %d not equal to expected %d\n", n, ns);
      					exit(-1);
    					} 	

    			fread(&Nchunk, size_int, 1, fid);


    			if(Nchunk>orgNvolume){
      					fprintf(stderr, "Error. ns = %d, nchunk = %d\n", ns, Nchunk);
      					fprintf(stderr, "System matrix chunk length %d is longer than Nvolume %d\n", Nchunk, orgNvolume);
      					exit(-1);
    					}

    			fread(Itemp,  size_int,   Nchunk, fid);
    			fread(SMtemp, size_float, Nchunk, fid);
	
	///////////////////main matrix multiplication/////////////////////////////////////////////////////////////
    		
   			for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++){
	   			if(singleProjFlag[rotationIndex][currentProjection]){
		    			if(fabs(sinogram[ns+rotationInSinogram[rotationIndex]])>1.0e-14){
			    			for(n=0;n<Nchunk;n++){
				   		 	for(basisIndex=0;basisIndex<basisNum;basisIndex++){
					    			coeffVolume[Itemp[n]+basisInVolume[basisIndex]] 
								+= sinogram[ns+rotationInSinogram[rotationIndex]] 
								* basis[basisIndex][singleProjInBasis[rotationIndex][currentProjection]] 
								* SMtemp[n];
							}}}}}

    				placeInCurrentProjection++;
    				if(placeInCurrentProjection==singleProjectionSize){
	   				placeInCurrentProjection=0;
	   				currentProjection++;
    					}

  		}// ns for loop closes here

  			for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++){
	  			free(singleProjInBasis[rotationIndex]);
	  			free(singleProjFlag[rotationIndex]);
  				}
 	 	
		free(singleProjInBasis);
  		free(singleProjFlag);
  		free(Itemp);
  		free(SMtemp);
  		free(basisInVolume);
 	 	fclose(fid);

}


///////////// Forward projection using a system-matrix, modified by a set of time basis functions ////////////
void GateBasisFP(int Nsinogram, int Nvolume, int singleProjectionSize, int basisNum, 
			int basisLen, int gatebasisNum, int gatebasisLen, float *sinogram, float *coeffVolume, float **basis, float **gatebasis,
			char *sysmat_fname){



  		int orgNvolume;
  		int ns, n, Nchunk, basisIndex, rotationIndex, rotationNum;
  		int currentProjection, placeInCurrentProjection, orgNsinogram, singleProjInRotation;
  		float *SMtemp;
  		int *Itemp,*basisInVolume,*rotationInSinogram;
  		int **singleProjFlag, **singleProjInBasis;
  		FILE *fid;
  		size_t size_int,size_intptr, size_float;

  		size_int = sizeof(int);
  		size_intptr = sizeof(int *);
  		size_float = sizeof(float);


 		 if(singleProjectionSize*basisLen!=Nsinogram)
  		{
    		fprintf(stderr, "Number of single projections %d and time frames %d does not match output sinogram size %d\n", 							singleProjectionSize,basisLen,Nsinogram);
    			exit(-1);
  		}

  		if((fid = fopen(sysmat_fname, "rb"))==NULL)
  		{
   			 fprintf(stderr, "Could not open sysmat file %s\n", sysmat_fname);
    			exit(-1);
  		}


  		fread(&orgNsinogram, size_int, 1, fid);

  		singleProjInRotation=0;
  		while(singleProjInRotation*singleProjectionSize<orgNsinogram)
	  		singleProjInRotation++;

  			currentProjection=0;
  			rotationNum=0;
  		while(currentProjection<basisLen)
  		{
	  		currentProjection+=singleProjInRotation;
	  		rotationNum++;
  		}
   		if(rotationNum>=MAX_ROTATIONS)
  		{
    			fprintf(stderr, "Using original sinogram size %d (from system matrix) for basis length %d exceeds %d rotatioms.\n", 				orgNsinogram, basisLen, MAX_ROTATIONS);
    			exit(-1);
  		}


  		fread(&orgNvolume, size_int, 1, fid);
  		if(Nvolume!=orgNvolume*basisNum)
  		{
    			fprintf(stderr, "Given Nvolume %d is not equal to expected %d * %d = %d.\n", 
				Nvolume, orgNvolume, basisNum, orgNvolume*basisNum);
    			exit(-1);
 		 }


  		Itemp = (int *)   malloc(size_int  *  orgNvolume);
  		SMtemp = (float *) malloc(size_float * orgNvolume);

  		basisInVolume = (int *)   malloc(size_int  *  basisNum);
  		for(basisIndex=0;basisIndex<basisNum;basisIndex++)
	  		basisInVolume[basisIndex]=basisIndex*orgNvolume;
  			rotationInSinogram = (int *)   malloc(size_int  *  rotationNum);
  		for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++)
	  		rotationInSinogram[rotationIndex]=rotationIndex*orgNsinogram;


  		singleProjInBasis = (int **)   malloc(size_intptr  *  rotationNum);
  		singleProjFlag = (int **)   malloc(size_intptr  *  rotationNum);

  		for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++)
  		{
	  		 singleProjInBasis[rotationIndex] = (int *)   malloc(size_int  *  singleProjInRotation);
	 		 singleProjFlag[rotationIndex] = (int *)   malloc(size_int  *  singleProjInRotation);

	  	for(n=0;n<singleProjInRotation;n++)
	  		{
				singleProjInBasis[rotationIndex][n] = rotationIndex*singleProjInRotation + n;

				if(singleProjInBasis[rotationIndex][n]<basisLen)
					singleProjFlag[rotationIndex][n] = 1;
				else
					singleProjFlag[rotationIndex][n] = 0;
	  		}
  		}

  		currentProjection=0;
  		placeInCurrentProjection=0;


  		for(ns=0;ns<orgNsinogram;ns++)
  		{

    			fread(&n, size_int, 1, fid);
   			 if(n!=ns)
    			{
      			fprintf(stderr, "Read in sinogram index %d not equal to expected %d\n", n, ns);
      			exit(-1);
    			}

    		fread(&Nchunk, size_int, 1, fid);
    		if(Nchunk>orgNvolume)
   		 {
      			fprintf(stderr, "Error. ns = %d, nchunk = %d\n", ns, Nchunk);
     			 fprintf(stderr, "System matrix chunk length %d is longer than Nvolume %d\n", Nchunk, orgNvolume);
      			exit(-1);
    		}

    		fread(Itemp,  size_int,   Nchunk, fid);
    		fread(SMtemp, size_float, Nchunk, fid);


    		for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++)
	    		if(singleProjFlag[rotationIndex][currentProjection])
		    		sinogram[ns+rotationInSinogram[rotationIndex]] = 0.;


    		for(n=0;n<Nchunk;n++)//for all voxel
	    		for(basisIndex=0;basisIndex<basisNum;basisIndex++)//for each basis function
    	   			 {
		   		if(fabs(coeffVolume[Itemp[n]+basisInVolume[basisIndex]])>1.0e-14)
			    		for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++)
		    	    			{

				    	if(singleProjFlag[rotationIndex][currentProjection])
					   sinogram[ns+rotationInSinogram[rotationIndex]] += 
						coeffVolume[Itemp[n]+basisInVolume[basisIndex]] * 
						basis[basisIndex][singleProjInBasis[rotationIndex][currentProjection]] * SMtemp[n];
			    }

	   	 }

    		placeInCurrentProjection++;
    		if(placeInCurrentProjection==singleProjectionSize)
    			{
	  		 placeInCurrentProjection=0;
	  		 currentProjection++;
   			 }

 		 }


  	for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++)
  		{
	  		free(singleProjInBasis[rotationIndex]);
	  		free(singleProjFlag[rotationIndex]);
  		}

  	free(singleProjInBasis);
  	free(singleProjFlag);
  	free(Itemp);
 	 free(SMtemp);
  	free(basisInVolume);

  	fclose(fid);

}


//////////////////////////////////////////////////////////////////////////////////////////////////
void TB_MLEM_recon(int Nsinogram, int Nvolume, int singleProjectionSize, int basisNum, int basisLen, 
		float *sinogram, float *coeffVolume, float **basis, char *sysmat_fname, 
		int LastIteration, int OutputNumber, int *IterationsForOutput, char *image_fname){
  
		FILE *fid;
  		int n, niter,index;
  		float *tmp_volume, *tmp_sinogram, *constant_denominator;
  		char temp_fname[128];
  		// --- initial volume assignment: all pixels are one
  		for(n=0; n<Nvolume; n++) coeffVolume[n] = 1.0;


		// ------ create ml_em variables:
  		tmp_sinogram = (float *) malloc((size_t) (sizeof(float) * Nsinogram));
  		constant_denominator = (float *) malloc((size_t) (sizeof(float) * Nvolume));
  		tmp_volume = (float *) malloc((size_t) (sizeof(float) * Nvolume));



		// --- compute  element-by-element inverse of efficiency matrix
  		for(n=0; n<Nsinogram; n++) tmp_sinogram[n] = 1.0;

  		TimeBasisBP(Nsinogram, Nvolume, singleProjectionSize, basisNum, basisLen, 
				tmp_sinogram, constant_denominator, basis, sysmat_fname);


  		for(n=0; n<Nvolume; n++)
    		if(constant_denominator[n] > EPSILON)  constant_denominator[n] = 1./ constant_denominator[n];
    		else constant_denominator[n] = ONE_OVER_EPSILON;


	//  -------- ITERATION LOOP --------
  		for(niter=1; niter<=LastIteration; niter++)
 			 {

    			//fprintf(stderr, "Iteration No %d of %d\n", niter, LastIteration);
			// compute the reprojection through the n-1 version of the file into tmp_sinogram
    		TimeBasisFP(Nsinogram, Nvolume, singleProjectionSize, basisNum,
				 basisLen, tmp_sinogram, coeffVolume, basis, sysmat_fname);


	// divide the sinogram by the tmp_sinogram
    	for(n=0; n<Nsinogram; n++)
      		if(sinogram[n] == 0.)  tmp_sinogram[n] = 0.;
      		else if(sinogram[n] < 0.)
        		/*nrerror("sinogram in MLEM smaller than zero");   */;
      		else
           	if(tmp_sinogram[n] > EPSILON) tmp_sinogram[n] = sinogram[n] / tmp_sinogram[n];
           	else                          tmp_sinogram[n] = sinogram[n] * ONE_OVER_EPSILON;
		// backproject the result into tmp_volume
      		TimeBasisBP(Nsinogram, Nvolume, singleProjectionSize, 
				basisNum, basisLen, tmp_sinogram, tmp_volume, basis, sysmat_fname);

		// multiply by the constant denominator
   		 for(n=0; n<Nvolume; n++)
    			{
      			coeffVolume[n] *=  constant_denominator[n] * tmp_volume[n];
      			if(coeffVolume[n] < 0.) coeffVolume[n] = 0.;
    			}



    		for(index=0;index<OutputNumber;index++)
			if(niter==IterationsForOutput[index])
    			{
				sprintf(temp_fname,"%s.iter%03d",image_fname,niter);
				if((fid = fopen(temp_fname, "wb"))==NULL)  //sprintf("%s.iter%d",image_fname,niter)
					{fprintf(stderr, "Could not open image file \"%s\" for writing\n", temp_fname);
					exit(1);}
			fwrite(coeffVolume, sizeof(float), Nvolume, fid);
			fclose(fid);
			}
  	}
	//end iteration here

  	// end: free up memory
 	 free(tmp_sinogram);
 	 free(tmp_volume);
  	 free(constant_denominator);

}



/*//////////////// Backward projection using a system-matrix, modified by a set of time basis functions *//////////////////
void TimeBasisBP(int Nsinogram, int Nvolume, int singleProjectionSize, 
			int basisNum, int basisLen, float *sinogram, float *coeffVolume, 
			float **basis, char *sysmat_fname){

  		int orgNvolume;
  		int ns, n, Nchunk, basisIndex, rotationIndex, rotationNum;
  		int currentProjection, placeInCurrentProjection, orgNsinogram, singleProjInRotation;
  		float *SMtemp;
  		int *Itemp,*basisInVolume,*rotationInSinogram;
  		int **singleProjFlag, **singleProjInBasis;
  		FILE *fid;
  		size_t size_int,size_intptr, size_float;

  		size_int = sizeof(int);
  		size_intptr = sizeof(int *);
  		size_float = sizeof(float);


  		/* Check if paramters match for the total sinogram size */
  		if(singleProjectionSize*basisLen!=Nsinogram)
  			{
    			fprintf(stderr, "Number of single projections %d and time frames %d does not match output sinogram size %d\n", 							singleProjectionSize,basisLen,Nsinogram);
    			exit(-1);
  			}
		/*open system matrix file*/
  		if((fid = fopen(sysmat_fname, "rb"))==NULL)
  			{
    				fprintf(stderr, "Could not open sysmat file %s\n", sysmat_fname);
    				exit(-1);
  			}


  	/* Read in static sinogram size from original system matrix */
  			fread(&orgNsinogram, size_int, 1, fid);

  			singleProjInRotation=0;
  		while(singleProjInRotation*singleProjectionSize<orgNsinogram)
	  		singleProjInRotation++;

  			currentProjection=0;
 			rotationNum=0;
  		while(currentProjection<basisLen)
  		{
	  		currentProjection+=singleProjInRotation;
	  		rotationNum++;
  		}
   		if(rotationNum>=MAX_ROTATIONS)
  		{
    		fprintf(stderr, "Using original sinogram size %d (from system matrix) for basis length %d exceeds %d rotatioms.\n", 
orgNsinogram, basisLen, MAX_ROTATIONS);
    		exit(-1);
  		}


  		/* Read in volume size, compare to expected */
  		fread(&orgNvolume, size_int, 1, fid);
  		if(Nvolume!=orgNvolume*basisNum)
 		 {
    			fprintf(stderr, "Given Nvolume %d is not equal to expected %d * %d = %d.\n", Nvolume, orgNvolume, basisNum, orgNvolume*basisNum);
    			exit(-1);
  		 }

 		 /* Set volume values to zero */
  		memset(coeffVolume, 0, size_float*Nvolume);
  		/* Allocate memory and compute variables to speed up calculations*/
  		Itemp = (int *)   malloc(size_int  *  orgNvolume);
 		SMtemp = (float *) malloc(size_float * orgNvolume);

  		basisInVolume = (int *)   malloc(size_int  *  basisNum);
  		for(basisIndex=0;basisIndex<basisNum;basisIndex++)
	  		basisInVolume[basisIndex]=basisIndex*orgNvolume;
  		rotationInSinogram = (int *)   malloc(size_int  *  rotationNum);
  		for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++)
	  		rotationInSinogram[rotationIndex]=rotationIndex*orgNsinogram;

  			singleProjInBasis = (int **)   malloc(size_intptr  *  rotationNum);
  			singleProjFlag = (int **)   malloc(size_intptr  *  rotationNum);

  		for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++)
  		{
	  		singleProjInBasis[rotationIndex] = (int *)   malloc(size_int  *  singleProjInRotation);
	  		singleProjFlag[rotationIndex] = (int *)   malloc(size_int  *  singleProjInRotation);

	  	for(n=0;n<singleProjInRotation;n++)
	  		{
				singleProjInBasis[rotationIndex][n] = rotationIndex*singleProjInRotation + n;
				if(singleProjInBasis[rotationIndex][n]<basisLen)
					singleProjFlag[rotationIndex][n] = 1;
				else
				singleProjFlag[rotationIndex][n] = 0;
	  		}
 		 }

  		currentProjection=0;
  		placeInCurrentProjection=0;

  		/* Main loop */
  		for(ns=0;ns<orgNsinogram;ns++)
  		{
    		/* Read volume index and compare to expected */
   		 fread(&n, size_int, 1, fid);
   		 if(n!=ns)
    		{
      			fprintf(stderr, "Read in sinogram index %d not equal to expected %d\n", n, ns);
      			exit(-1);
    		}


   	 /* Fread chunk size and indices and SM chunks */
    		fread(&Nchunk, size_int, 1, fid);
    		if(Nchunk>orgNvolume)
    		{
      			fprintf(stderr, "Error. ns = %d, nchunk = %d\n", ns, Nchunk);
      			fprintf(stderr, "System matrix chunk length %d is longer than Nvolume %d\n", Nchunk, orgNvolume);
      			exit(-1);
    		}
    		fread(Itemp,  size_int,   Nchunk, fid);
    		fread(SMtemp, size_float, Nchunk, fid);

    		/* Do loop multiplication, going through volume voxels, basis functions, and rotations */
   		 for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++)
	   		/*If this time point exists, compute the value for the sinogram */
    	    	if(singleProjFlag[rotationIndex][currentProjection])
		    if(fabs(sinogram[ns+rotationInSinogram[rotationIndex]])>1.0e-14)
			    for(n=0;n<Nchunk;n++)
				    for(basisIndex=0;basisIndex<basisNum;basisIndex++)
					    coeffVolume[Itemp[n]+basisInVolume[basisIndex]] 
					+= sinogram[ns+rotationInSinogram[rotationIndex]] 
					* basis[basisIndex][singleProjInBasis[rotationIndex][currentProjection]] * SMtemp[n];

    		placeInCurrentProjection++;
    		if(placeInCurrentProjection==singleProjectionSize)
    			{
	   			placeInCurrentProjection=0;
	   			currentProjection++;
    			}

  		}

  	/* Free memory, close file and return */
  		for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++)
  		{
	  		free(singleProjInBasis[rotationIndex]);
	  		free(singleProjFlag[rotationIndex]);
  		}

		
  	free(singleProjInBasis);
  	free(singleProjFlag);
  	free(Itemp);
  	free(SMtemp);
  	free(basisInVolume);
 	 fclose(fid);


}


/*///////////// Forward projection using a system-matrix, modified by a set of time basis functions */////////////
void TimeBasisFP(int Nsinogram, int Nvolume, int singleProjectionSize, int basisNum, 
			int basisLen, float *sinogram, float *coeffVolume, float **basis, 
				char *sysmat_fname){

  		int orgNvolume;
  		int ns, n, Nchunk, basisIndex, rotationIndex, rotationNum;
  		int currentProjection, placeInCurrentProjection, orgNsinogram, singleProjInRotation;
  		float *SMtemp;
  		int *Itemp,*basisInVolume,*rotationInSinogram;
  		int **singleProjFlag, **singleProjInBasis;
  		FILE *fid;
  		size_t size_int,size_intptr, size_float;

  		size_int = sizeof(int);
  		size_intptr = sizeof(int *);
  		size_float = sizeof(float);

  		/* Check if paramters match for the total sinogram size */
 		 if(singleProjectionSize*basisLen!=Nsinogram)
  		{
    		fprintf(stderr, "Number of single projections %d and time frames %d does not match output sinogram size %d\n", 							singleProjectionSize,basisLen,Nsinogram);
    			exit(-1);
  		}

  		if((fid = fopen(sysmat_fname, "rb"))==NULL)
  		{
   			 fprintf(stderr, "Could not open sysmat file %s\n", sysmat_fname);
    			exit(-1);
  		}

 		 /* Read in static sinogram size from original system matrix */
  		fread(&orgNsinogram, size_int, 1, fid);

  		singleProjInRotation=0;
  		while(singleProjInRotation*singleProjectionSize<orgNsinogram)
	  		singleProjInRotation++;

  			currentProjection=0;
  			rotationNum=0;
  		while(currentProjection<basisLen)
  		{
	  		currentProjection+=singleProjInRotation;
	  		rotationNum++;
  		}
   		if(rotationNum>=MAX_ROTATIONS)
  		{
    			fprintf(stderr, "Using original sinogram size %d (from system matrix) for basis length %d exceeds %d rotatioms.\n", 				orgNsinogram, basisLen, MAX_ROTATIONS);
    			exit(-1);
  		}


 		 /* Read in volume size, compare to expected */
  		fread(&orgNvolume, size_int, 1, fid);
  		if(Nvolume!=orgNvolume*basisNum)
  		{
    			fprintf(stderr, "Given Nvolume %d is not equal to expected %d * %d = %d.\n", 
				Nvolume, orgNvolume, basisNum, orgNvolume*basisNum);
    			exit(-1);
 		 }

  		/* Allocate memory and compute variables to speed up calculations*/
  		Itemp = (int *)   malloc(size_int  *  orgNvolume);
  		SMtemp = (float *) malloc(size_float * orgNvolume);

  		basisInVolume = (int *)   malloc(size_int  *  basisNum);
  		for(basisIndex=0;basisIndex<basisNum;basisIndex++)
	  		basisInVolume[basisIndex]=basisIndex*orgNvolume;
  			rotationInSinogram = (int *)   malloc(size_int  *  rotationNum);
  		for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++)
	  		rotationInSinogram[rotationIndex]=rotationIndex*orgNsinogram;


  		singleProjInBasis = (int **)   malloc(size_intptr  *  rotationNum);
  		singleProjFlag = (int **)   malloc(size_intptr  *  rotationNum);

  		for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++)
  		{
	  		 singleProjInBasis[rotationIndex] = (int *)   malloc(size_int  *  singleProjInRotation);
	 		 singleProjFlag[rotationIndex] = (int *)   malloc(size_int  *  singleProjInRotation);

	  	for(n=0;n<singleProjInRotation;n++)
	  		{
				singleProjInBasis[rotationIndex][n] = rotationIndex*singleProjInRotation + n;

				if(singleProjInBasis[rotationIndex][n]<basisLen)
					singleProjFlag[rotationIndex][n] = 1;
				else
					singleProjFlag[rotationIndex][n] = 0;
	  		}
  		}

  		currentProjection=0;
  		placeInCurrentProjection=0;

  		/* Main loop */
  		for(ns=0;ns<orgNsinogram;ns++)
  		{
    		/* Read volume index and compare to expected */
    			fread(&n, size_int, 1, fid);
   			 if(n!=ns)
    			{
      			fprintf(stderr, "Read in sinogram index %d not equal to expected %d\n", n, ns);
      			exit(-1);
    			}
    		/* Fread chunk size and indices and SM chunks */
    		fread(&Nchunk, size_int, 1, fid);
    		if(Nchunk>orgNvolume)
   		 {
      			fprintf(stderr, "Error. ns = %d, nchunk = %d\n", ns, Nchunk);
     			 fprintf(stderr, "System matrix chunk length %d is longer than Nvolume %d\n", Nchunk, orgNvolume);
      			exit(-1);
    		}

    		fread(Itemp,  size_int,   Nchunk, fid);
    		fread(SMtemp, size_float, Nchunk, fid);

    		/* Set the correct sinogram values to zero */
    		for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++)
	    		if(singleProjFlag[rotationIndex][currentProjection])
		    		sinogram[ns+rotationInSinogram[rotationIndex]] = 0.;


    		/* Do loop multiplication, going through volume voxels, basis functions, and rotations */
    		for(n=0;n<Nchunk;n++)//for all voxel
	    		for(basisIndex=0;basisIndex<basisNum;basisIndex++)//for each basis function
    	   			 {
		   		if(fabs(coeffVolume[Itemp[n]+basisInVolume[basisIndex]])>1.0e-14)
			    		for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++)
		    	    			{
				   	 /* If this time point exists, compute the value for the sinogram */
				    	if(singleProjFlag[rotationIndex][currentProjection])
					   sinogram[ns+rotationInSinogram[rotationIndex]] += 
						coeffVolume[Itemp[n]+basisInVolume[basisIndex]] * 
						basis[basisIndex][singleProjInBasis[rotationIndex][currentProjection]] * SMtemp[n];
			    }

	   	 }

    		placeInCurrentProjection++;
    		if(placeInCurrentProjection==singleProjectionSize)
    			{
	  		 placeInCurrentProjection=0;
	  		 currentProjection++;
   			 }

 		 }

 	 /* Free memory, close file and return */
  	for(rotationIndex=0;rotationIndex<rotationNum;rotationIndex++)
  		{
	  		free(singleProjInBasis[rotationIndex]);
	  		free(singleProjFlag[rotationIndex]);
  		}

  	free(singleProjInBasis);
  	free(singleProjFlag);
  	free(Itemp);
 	 free(SMtemp);
  	free(basisInVolume);

  	fclose(fid);

}



