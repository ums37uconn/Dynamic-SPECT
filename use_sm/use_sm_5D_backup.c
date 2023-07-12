#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


// ----- This block defines implementation-specific parameters -----
// 	This define should be on: for now I only have ray-driven system matrices
#define RAY_DRIVEN

// input sinogram can be either 16-bit unsigned integer (unsigned short) or float
#define FLOAT_INPUT_SINOGRAM
// uncomment this, if input sinogram is float, leave like this otherwise

// Volume is best output as 4-bit float
#define FLOAT_INPUT_VOLUME

// No of iterations
// #define N_ITERATIONS 25

// Max iterations and outputs of different iterations
#define MAX_ITERATIONS 1000
#define MAX_OUTPUTS 100
// Max number of rotations of camera accepted by the forward/backward projection with a time basis
#define	MAX_ROTATIONS	1000

// This defines how small zero is, and what is one divided by zero
#define ONE_OVER_EPSILON 1.0e+13
#define EPSILON 1.0e-13

// --- functions.  
//This contains everything we need except for nrutil functions

	int main(int argc, char **argv);
	void MLEM_recon(int Nsinogram, float *sinogram, int Nvolume, 
		float *volume, char *fname, int LastIteration, int OutputNumber, 
		int *IterationsForOutput, char *image_fname);
	void VolDrBP_SM(int Nsinogram, float *sinogram, int Nvolume, float *volume, char *fname);
	void VolDrFP_SM(int Nsinogram, float *sinogram, int Nvolume, float *volume, char *fname);
	void RayDrBP_SM(int Nsinogram, float *sinogram, int Nvolume, float *volume, char *fname);
	void RayDrFP_SM(int Nsinogram, float *sinogram, int Nvolume, float *volume, char *fname);

//.......for dynamic recon

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

//........for gate + dynamic recon

	void GB_MLEM_recon(int Nsinogram, int Nvolume, int singleProjectionSize, 
			int basisNum, int basisLen, float *sinogram, float *coeffVolume, 
			float **basis, char *sysmat_fname, int LastIteration, int OutputNumber, int *IterationsForOutput, char *image_fname);
	void GateBasisBP(int Nsinogram, int Nvolume, int singleProjectionSize, int basisNum, 
			int basisLen, float *sinogram, float *coeffVolume, float **basis, 
			char *sysmat_fname);
	void GateBasisFP(int Nsinogram, int Nvolume, int singleProjectionSize, 
			int basisNum, int basisLen, float *sinogram, float *coeffVolume, 
			float **basis, char *sysmat_fname);

// ------------ END HEADER, BEGIN FUNCTIONS ---------

// --- MAIN CODE.
//     command line input arguments are explained in error message printed when a.out is started
	
	int main(int argc, char **argv){
  
		int n, Nvolume, Nsinogram;
  		int IterationsForOutput[MAX_OUTPUTS];
  		int OutputNumber=1;
  		int index;

 		int tbFlag,tbNum,tbLen,singleProjectionSize;
		int gbFlag,gbNum,gbLen;
  		
		FILE *fid;
  		float **basis;
  		float *sinogram, *volume;

  		char sino_fname[128], image_fname[128], sysmat_fname[128], basis_fname[128], Action;
  		unsigned short tmp;

// DEBUG: on
// Fares Alhassen, UCSF
// adding variables to replace magic number used to control input parameters

  		int minNumberOfArguments = 6;
 		int maxNumberOfArguments = 9;

//Error messages

  		if((argc < minNumberOfArguments ) || (argc > maxNumberOfArguments + MAX_OUTPUTS)){
   			 fprintf(stderr, "\nUsage:\n%s sysmat.file sinogram.file image.file" 
				"operation n_iterations {timebasis.file SingleProjectionSize}"
 				"(Iteration numbers for output)\n",argv[0]);
			
    			 fprintf(stderr, "Where \"operation\" can be either R for reconstruction," 
				"or F or B for forward or backprojection (CAPITAL LETTERS).\n");
    			 fprintf(stderr, "Or r\\f\\b (non capital letters) for these operations"
				 "with a time basis file {following should be a time basis file," 
				"and the number of data points or projection from a single angle}.\n");
    			 fprintf(stderr, "(Iteration numbers can be one or more number of iterations" 
				"below %d, such as: 10 25 100.)\n",MAX_ITERATIONS);
   		 	exit(1);
 		 }

 		strcpy(sysmat_fname, argv[1]);
  		strcpy(sino_fname,   argv[2]);
 		strcpy(image_fname,  argv[3]);
 		Action = argv[4][0];

  		switch(Action){
	
    			case 'g': Action = 'g'; tbFlag = 1; gbFlag=1; break;//reconstruction with gate and time basis on			
			
			case 'r': Action = 'r'; tbFlag = 1; break;//reconstruction with time basis on
    			case 'R': Action = 'r'; tbFlag = 0; break;//reconstruction time basis off
   			case 'f': Action = 'f'; tbFlag = 1; break;//forward projection  with time  basis on
   			case 'F': Action = 'f'; tbFlag = 0; break;//forward projection with time basis off
    			case 'b': Action = 'b'; tbFlag = 1; break;//backward projection with time basis on
    			case 'B': Action = 'b'; tbFlag = 0; break;//backward projection with time basis off
   			default:
      			fprintf(stderr, "Action setting should be: \"R\" or \"r\"  for reconstruction, \"F\" or \"f\" for forwardprojection,\n");
      			fprintf(stderr, "or \"B\" or \"b\" for backprojection, unrecongnized option \"%s\"\n", argv[4]);
     			exit(1);

 			 }

 	 // get number of iterations to compute
 		int LastIteration = atoi( argv[ minNumberOfArguments - 1 ] );

  	// -- read sinogram and image size from the sysmat file.
  	// DEBUG: off

//read system matrix file
  	// fprintf( stderr, "Reading sinogram and image size from system matrix file...\n" );
 	 	if (( fid = fopen(sysmat_fname, "rb" )) == NULL ) {
      			fprintf(stderr, "Could not open sysmat file \"%s\" for reading\n", sysmat_fname);
      			exit(1);
  		}


  		fread( &Nsinogram, sizeof( int ), 1, fid );
  		fread( &Nvolume, sizeof( int ), 1, fid );

  		fclose( fid );

/////////////////////////////////////////////////////////////////////////////////////////////////////////

//If using time basis file, than read it //for 4D

  	if ( tbFlag ){ 

	  	if (argc > 6){
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
	 		 }

	  	else{
			fprintf(stderr, "After operation r\\f\\b should come the time basis file," 
				"and following the number of data points of projection from a single angle \n");
			exit(1);
	  	    }

	 		Nsinogram = singleProjectionSize * tbLen;
	  		Nvolume = Nvolume * tbNum;

 	      }
	//tbFlag closes here


////////////////////////////////////////////////////////////////////////////////////////////

  		// -- check if we are given the number of iterations
  	
		if(argc > minNumberOfArguments + tbFlag * 2 ){
			LastIteration=0;
			OutputNumber = argc - minNumberOfArguments - tbFlag * 2;
			for(index=0; index<OutputNumber; index++){
				IterationsForOutput[index]=atoi(argv[index+minNumberOfArguments+tbFlag*2]);
				if((IterationsForOutput[index]<=0) || (IterationsForOutput[index]>MAX_ITERATIONS)){
					fprintf(stderr, "Iteration number %d for output must be an integer below %d.\n",
						atoi(argv[index+minNumberOfArguments+tbFlag*2]),MAX_ITERATIONS);
   					exit(1);
					}
			if(IterationsForOutput[index]>LastIteration)
					LastIteration=IterationsForOutput[index];
			}}

  		else{
	  		IterationsForOutput[0]=LastIteration;
		    }


 	 	// allocations
 	 	// DEBUG: off
  	 	// fprintf( stderr, "Allocating sinogram with %d elements...\n", Nsinogram );
  
			sinogram = (float *) malloc(sizeof(float)*Nsinogram);

  		// fprintf( stderr, "Allocating volume with %d elements...\n", Nvolume );

  			volume   = (float *) malloc(sizeof(float)*Nvolume);


 		// read input data
  		if(Action=='f')  // read volume
  		{
    				if((fid = fopen(image_fname, "rb"))==NULL){
      				fprintf(stderr, "Could not open image data file \"%s\" for reading\n", image_fname);
      				exit(1);
   						 }

			#ifndef FLOAT_INPUT_VOLUME
    			for(n=0; n<Nvolume; n++){ 
				fread(&tmp, 2, 1, fid); // 2 = sizeof(unsigned short) == 16 bit
      				volume[n] = tmp;
					}
			#else
   				 fread(volume, sizeof(float), Nvolume, fid);
			#endif
    				fclose(fid);

 		 }

 		else //-----------read sinogram
 		{
    				// DEBUG: off
   		 		// fprintf( stderr, "Reading sinogram with %d elements...\n", Nsinogram );
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

  		}

////////////////////////////////////////////////////////////////////////////////////////////////////

// At this point, we have read all the required input data, except sysmat

  	//if(tbFlag)fprintf(stderr, "Using %d basis functions of length %d.\n", tbNum, tbLen);

  	//  --- Main action ---
  		switch(Action)
  			{
    				case 'r':
   	//fprintf(stderr, "Reconstruction volume of length %d from sinogram of length %d\n", Nvolume, Nsinogram);
      				if(tbFlag)
	      				TB_MLEM_recon(Nsinogram, Nvolume, singleProjectionSize, 
					tbNum, tbLen, sinogram, volume, basis, sysmat_fname, 
					LastIteration, OutputNumber, IterationsForOutput, image_fname);
      				else
	      				MLEM_recon(Nsinogram, sinogram, Nvolume, volume, sysmat_fname, 
					LastIteration, OutputNumber, IterationsForOutput, image_fname);
      					break;
      
    				case 'f':
      				fprintf(stderr, "Forward project volume of length %d into sinogram of length %d\n", 
					Nvolume, Nsinogram);
     		 			if(tbFlag)
	    					TimeBasisFP(Nsinogram, Nvolume, singleProjectionSize, tbNum, tbLen, 
						sinogram, volume, basis, sysmat_fname);
      						else

					#ifdef RAY_DRIVEN
      						RayDrFP_SM(Nsinogram, sinogram, Nvolume, volume, sysmat_fname);
					#else
      						VolDrFP_SM(Nsinogram, sinogram, Nvolume, volume, sysmat_fname);
					#endif
     					 break;

    				case 'b':
     				fprintf(stderr, "Backproject volume of length %d from sinogram of length %d\n",
					 Nvolume, Nsinogram);
      				if(tbFlag)
	    			TimeBasisBP(Nsinogram, Nvolume, singleProjectionSize, tbNum, tbLen, sinogram, 
					volume, basis, sysmat_fname);
      				else

					#ifdef RAY_DRIVEN
     	 					RayDrBP_SM(Nsinogram, sinogram, Nvolume, volume, sysmat_fname);
					#else
     					 	VolDrBP_SM(Nsinogram, sinogram, Nvolume, volume, sysmat_fname);
					#endif
      					break;
    					default:
      				/*nrerror("This point should not be reachable in use_sm!");*/  ;
  		}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// --- Save results ---
  		if(Action=='f') // save sinogram
  			{
    				if((fid = fopen(sino_fname, "wb"))==NULL)
    					{ 
					fprintf(stderr, "Could not open sinogram file \"%s\" for writing\n", sino_fname);
      					exit(1);
					}
    				fwrite(sinogram, sizeof(float), Nsinogram, fid);
    			fclose(fid);
  			}

  		else if(Action=='b')  // save volume
 	 		{
   			 if((fid = fopen(image_fname, "wb"))==NULL){ 
				fprintf(stderr, "Could not open image file \"%s\" for writing\n", image_fname);
     		 		exit(1);
				}
    			fwrite(volume, sizeof(float), Nvolume, fid);
    			fclose(fid);
  			}

	
  	//free memory and return
  			free(volume);
  			free(sinogram);

  
  	if(tbFlag){
	 	for(index=0;index<tbNum;index++)
			 free(basis[index]);
			 free(basis);
 		 }
 	 return 1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MLEM_recon(int Nsinogram, float *sinogram, int Nvolume, 
		float *volume, char *fname, int LastIteration, int OutputNumber, 
		int *IterationsForOutput, char *image_fname){

  		FILE *fid;
  		int n, niter,index;
  		float *tmp_volume, *tmp_sinogram, *constant_denominator;

		float *volume_old;

  		char temp_fname[128];
  		// --- initial volume assignment: all pixels are one
  		for(n=0; n<Nvolume; n++) volume[n] = 1.0;

		// ------ create ml_em variables:
  			tmp_sinogram = (float *) malloc((size_t) (sizeof(float) * Nsinogram));
  			constant_denominator = (float *) malloc((size_t) (sizeof(float) * Nvolume));
  			tmp_volume = (float *) malloc((size_t) (sizeof(float) * Nvolume));
  			volume_old = (float *) malloc((size_t) (sizeof(float) * Nvolume));
		// --- compute  element-by-element inverse of efficiency matrix
 			 for(n=0; n<Nsinogram; n++) tmp_sinogram[n] = 1.;
  
		#ifdef RAY_DRIVEN
      			RayDrBP_SM(Nsinogram, tmp_sinogram, Nvolume, constant_denominator, fname);
		#else
     			VolDrBP_SM(Nsinogram, tmp_sinogram, Nvolume, constant_denominator, fname);
		#endif

 		for(n=0; n<Nvolume; n++)
    			if(constant_denominator[n] > EPSILON)  
					constant_denominator[n] = 1./ constant_denominator[n];
    			else constant_denominator[n] = ONE_OVER_EPSILON;


	//  -------- ITERATION LOOP --------
  		for(niter=1; niter<=LastIteration; niter++){
	 			for(n=0; n<Nvolume; n++)
					volume_old[n]=volume[n];

			//fprintf(stderr, "Iteration No %d of %d\n", niter, LastIteration);
		//compute the reprojection through the n-1 version of the file into tmp_sinogram
		#ifdef RAY_DRIVEN
      			RayDrFP_SM(Nsinogram, tmp_sinogram, Nvolume,  volume, fname);
		#else
      			VolDrFP_SM(Nsinogram, tmp_sinogram, Nvolume,  volume, fname);
		#endif

		//divide the sinogram by the tmp_sinogram
    		for(n=0; n<Nsinogram; n++)
     		 	if(sinogram[n] == 0.)  tmp_sinogram[n] = 0.;
      			else if(sinogram[n] < 0.)
        		/*nrerror("sinogram in MLEM smaller than zero");   */;
      			else

           	if(tmp_sinogram[n] > EPSILON) tmp_sinogram[n] = sinogram[n] / tmp_sinogram[n];
          	else                          tmp_sinogram[n] = sinogram[n] * ONE_OVER_EPSILON;
           

		// backproject the result inot tmp_volume
		#ifdef RAY_DRIVEN
      			RayDrBP_SM(Nsinogram, tmp_sinogram, Nvolume,  tmp_volume, fname);
		#else
      			VolDrBP_SM(Nsinogram, tmp_sinogram, Nvolume,  tmp_volume, fname);
		#endif

		// multiply by the constant denominator
    		for(n=0; n<Nvolume; n++)
    			{
      			volume[n] *=  constant_denominator[n] * tmp_volume[n];
      			if(volume[n] < 0.0) volume[n] = 0.0;
    			}

		// large buffer size
    			int outputBuffersize = strlen( temp_fname ) + 1;
	//    if ( strlen( temp_fname ) + 1 > outputBuffersize ) {
	//                fprintf( stderr, "Could not open image file \"%s\" for writing: filename is too long for buffer.\n", temp_fname );
	//                exit( EXIT_FAILURE );
	//            }

		float temp=0.0;
 		for(n=0; n<Nvolume; n++)temp+=(volume[n]-volume_old[n])*(volume[n]-volume_old[n]);
	 		//fprintf(stderr, "iteration number=%d      chi square=%f\n", niter, temp);
			fprintf(stderr, "%f\n", temp);


    		for( index = 0; index < OutputNumber; index++ )
        			if( niter == IterationsForOutput[ index ] ) {
           	 // DEBUG: on
            	// edit FA / UCSF 2011.06.02
            	// causes buffer overflow error if input filepath too large
            		sprintf( temp_fname, "%s.iter%03d", image_fname, niter );
            		// snprintf( temp_fname, outputBuffersize, "%s.iter%03d", image_fname, niter );
            		if ( ( fid = fopen( temp_fname, "wb" ) ) == NULL )  //sprintf("%s.iter%d",image_fname,niter)
                		{ fprintf(stderr, "Could not open image file \"%s\" for writing\n", temp_fname);
               			 exit(1);}
           			 fwrite(volume, sizeof(float), Nvolume, fid);
            			fclose(fid);
        		}

  		}

  		// end: free up memory
 			 free(tmp_sinogram);
  			free(tmp_volume);
  			free(volume_old);
  			free(constant_denominator);
  			return;

}


/////////////// Voxel driven system-matrix based backprojection///////////////////////////////////////////
void VolDrBP_SM(int Nsinogram, float *sinogram, int Nvolume, float *volume, char *fname){
  
		int nv, n, Nchunk;
  		float *SMtemp;
  		int *Itemp;
  		FILE *fid;
  		size_t size_int, size_float;

  		size_int = sizeof(int);
  		size_float = sizeof(float);

  		if((fid = fopen(fname, "rb"))==NULL)
  		{
    			fprintf(stderr, "Could not open sysmat file %s\n", fname);
    			exit(1);
  		}
  	// read in Nsinogram, compare to given
  		fread(&n, size_int, 1, fid);
  			if(n!=Nsinogram)
  				{
    					fprintf(stderr, "Read in Nvsinogram %d not equal to expected %d\n", n, Nsinogram);
   					exit(1);
 				 }
  	// read in Nvolume, compare to given
  		fread(&n, size_int, 1, fid);
  		if(n!=Nvolume)
  			{
    				fprintf(stderr, "Read in Nvolume %d not equal to expected %d\n", n, Nvolume);
    				exit(1);
  			}
 	 // set volume values to zero
 	 // memset(volume, 0, size_float*Nvolume); // don't need to do it here, do it in the loop
  	// Declare sm chunk and volume index chunk
  		Itemp = (int *)   malloc(size_int    * Nsinogram);
  		SMtemp = (float *) malloc(size_float * Nsinogram);

  	// Main loop
  	for(nv=0;nv<Nvolume;nv++)
  	{
    		// read volume index and compare to expected
    		fread(&n, size_int, 1, fid);
    		if(n!=nv)
    		{
      			fprintf(stderr, "Read in voxel index %d not equal to expected %d\n", n, nv);
     			 exit(1);
    		}
    	// fread chunk size and indices and SM chunks
    		fread(&Nchunk, size_int, 1, fid);
    		if(Nchunk>Nsinogram)
    			{
      			fprintf(stderr, "Sinogram chunk %d is longer than Nsinogram %d\n", Nchunk, Nsinogram);
     			 exit(1);
    			}
    		fread(Itemp, size_int,    Nchunk, fid);
    		fread(SMtemp, size_float, Nchunk, fid);
    	// do loop muptiplication
    		volume[nv]=0.;
   		 for(n=0;n<Nchunk;n++)
     		 if(fabs(sinogram[Itemp[n]])>1.0e-11)
       			 volume[nv] += sinogram[Itemp[n]]*SMtemp[n];

 	 }

  	// done, free memory, close file and return
  	 fclose(fid);
  	 free(Itemp);
 	 free(SMtemp);
  	 return;
}


/////////////// Voxel driven system-matrix based forward projection///////////////////////////////////////////
void VolDrFP_SM(int Nsinogram, float *sinogram, int Nvolume, float *volume, char *fname){
  
		int nv, n, Nchunk;
  		float *SMtemp;
  		int *Itemp;
 		 FILE *fid;
  		size_t size_int, size_float;

  		size_int = sizeof(int);
 		 size_float = sizeof(float);


  		if((fid = fopen(fname, "rb"))==NULL)
  			{
    			fprintf(stderr, "Could not open sysmat file %s\n", fname);
    			exit(1);
 			 }
  		// read in Nsinogram, compare to given
 		 fread(&n, size_int, 1, fid);
  		if(n!=Nsinogram)
 			 {
   			 fprintf(stderr, "Read in Nvsinogram %d not equal to expected %d\n", n, Nsinogram);
   			 exit(1);
  			}
 		 // read in Nvolume, compare to given
  		fread(&n, size_int, 1, fid);
  		if(n!=Nvolume)
 		 {
    			fprintf(stderr, "Read in Nvolume %d not equal to expected %d\n", n, Nvolume);
   			 exit(1);
  		 }

  	// set sinogram values to zero
  		memset(sinogram, 0, size_float*Nsinogram);
  	// Declare sm chunk and volume index chunk
 	 	Itemp = (int *)   malloc(size_int  * Nsinogram);
  		SMtemp = (float *) malloc(size_float*Nsinogram);

  	// Main loop
  	for(nv=0;nv<Nvolume;nv++)
  		{
    		// read volume index and compare to expected
    			fread(&n, size_int, 1, fid);
    			if(n!=nv)
    			{
      			fprintf(stderr, "Read in voxel index %d not equal to expected %d\n", n, nv);
      			exit(1);
    			}
    		// fread chunk size and indices and SM chunks
    			fread(&Nchunk, size_int, 1, fid);
    			if(Nchunk>Nsinogram)
    			{
      			fprintf(stderr, "Sinogram chunk %d is longer than Nsinogram %d\n", Nchunk, Nsinogram);
      			exit(1);
    			}

    			fread(Itemp,  size_int,   Nchunk, fid);
    			fread(SMtemp, size_float, Nchunk, fid);
   		 // if voxel is not zero, do loop muptiplication
    			if(fabs(volume[nv])>1.0e-11)
      		for(n=0;n<Nchunk;n++)
        		sinogram[Itemp[n]] += volume[nv]*SMtemp[n];
  	}
  	// done, free memory, close file and return
  	fclose(fid);
  	free(Itemp);
  	free(SMtemp);
  	return;
}


/////////////// Sinogram driven system-matrix based backprojection////////////////////////////////
void RayDrBP_SM(int Nsinogram, float *sinogram, int Nvolume, float *volume, char *fname){
  
		int ns, n, Nchunk;
  		float *SMtemp;
  		int *Itemp;
  		FILE *fid;
  		size_t size_int, size_float;

  		size_int = sizeof(int);
  		size_float = sizeof(float);


  		if((fid = fopen(fname, "rb"))==NULL)
  			{
    				fprintf(stderr, "Could not open sysmat file %s\n", fname);
    				exit(1);
 			 }
  		// read in Nsinogram, compare to given
  		fread(&n, size_int, 1, fid);
  		if(n!=Nsinogram)
  			{
    				fprintf(stderr, "Read in Nvsinogram %d not equal to expected %d\n", n, Nsinogram);
   				 exit(1);
 			 }
 		// read in Nvolume, compare to given
  		fread(&n, size_int, 1, fid);
  		if(n!=Nvolume)
 		 	{
   			 fprintf(stderr, "Read in Nvolume %d not equal to expected %d\n", n, Nvolume);
    			 exit(1);
  			}

  		// set volume values to zero
  			memset(volume, 0, size_float*Nvolume);
  		// Declare sm chunk and volume index chunk
  			Itemp = (int *)   malloc(size_int  * Nvolume);
  			SMtemp = (float *) malloc(size_float*Nvolume);

  		// Main loop
  		for(ns=0;ns<Nsinogram;ns++)
  			{
    		// read volume index and compare to expected
    			fread(&n, size_int, 1, fid);
    			if(n!=ns)
    			{
      				fprintf(stderr, "Read in sinogram index %d not equal to expected %d\n", n, ns);
      				exit(1);
    			}
    		// fread chunk size and indices and SM chunks
    			fread(&Nchunk, size_int, 1, fid);
    			if(Nchunk>Nvolume)
    			{
      			fprintf(stderr, "Sinogram chunk %d is longer than Nvolume %d\n", Nchunk, Nvolume);
      			exit(1);
    			}

    		fread(Itemp, size_int,    Nchunk, fid);
    		fread(SMtemp, size_float, Nchunk, fid);
    	// is sinogram pixel is non-zero, do loop muptiplication
   		 if(fabs(sinogram[ns])>1.0e-14)
      		for(n=0;n<Nchunk;n++)
       			 volume[Itemp[n]] += sinogram[ns]*SMtemp[n];

 	}

  	// done, free memory, close file and return
 	fclose(fid);
  	free(Itemp);
  	free(SMtemp);
  	return;
}


/////////////// Sinogram driven system-matrix based forward projection////////////////////////////////////////
void RayDrFP_SM(int Nsinogram, float *sinogram, int Nvolume, float *volume, char *fname){
  
		int ns, n, Nchunk;
  		float *SMtemp;
  		int *Itemp;
 		 FILE *fid;
  		size_t size_int, size_float;

 		 size_int = sizeof(int);
  		size_float = sizeof(float);


  		if((fid = fopen(fname, "rb"))==NULL)
  		{
    			fprintf(stderr, "Could not open sysmat file %s\n", fname);
    			exit(1);
  		}
  		// read in Nsinogram, compare to given
  		fread(&n, size_int, 1, fid);
  		if(n!=Nsinogram)
  			{
    			fprintf(stderr, "Read in Nvsinogram %d not equal to expected %d\n", n, Nsinogram);
    			exit(1);
  			}
  		// read in Nvolume, compare to given
  		fread(&n, size_int, 1, fid);
  		if(n!=Nvolume)
  			{
    			fprintf(stderr, "Read in Nvolume %d not equal to expected %d\n", n, Nvolume);
    			exit(1);
  			}
  		// set sinogram values to zero
 		 // memset(sinogram, 0, size_float*Nsinogram); // don't need to do it here, done in the loop
 		 // Declare sm chunk and volume index chunk
  		Itemp = (int *)   malloc(size_int  *  Nvolume);
  		SMtemp = (float *) malloc(size_float* Nvolume);

  	// Main loop
  		for(ns=0;ns<Nsinogram;ns++)
  		{
    			// read volume index and compare to expected
   			 fread(&n, size_int, 1, fid);
    			if(n!=ns)
    			{
      			fprintf(stderr, "Read in sinogram index %d not equal to expected %d\n", n, ns);
      			exit(1);
    			}
   	 // fread chunk size and indices and SM chunks
    		fread(&Nchunk, size_int, 1, fid);
    		if(Nchunk>Nvolume)
    		{
      			fprintf(stderr, "ns = %d, nchunk = %d\n", ns, Nchunk);
      			fprintf(stderr, "Forward proj: System matrix chunk length %d is longer than Nvolume %d\n", Nchunk, Nvolume);
     			 exit(1);
    		}
    		fread(Itemp,  size_int,   Nchunk, fid);
    		fread(SMtemp, size_float, Nchunk, fid);
    		sinogram[ns] = 0.;
   	 // do loop muptiplication
    		for(n=0;n<Nchunk;n++)
      		if(fabs(volume[Itemp[n]])>1.0e-14)
        		sinogram[ns] += volume[Itemp[n]]*SMtemp[n];

  	}

  	// done, free memory, close file and return
 	 fclose(fid);
  	free(Itemp);
  	free(SMtemp);
  	return;
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



