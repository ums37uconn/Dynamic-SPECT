/* DEBUG: on */
/* comment off below statement for 32-bit system */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cfiMesh2.h"

#define SPARSE_FORMAT
#define PARM_NUM		24
#define MAX_STR_LEN		1024

cfiDebugStatus g_debugStatus = (cfiProgressDebug|cfiArgDebug|cfiReturnDebug);
/* cfiDebugStatus g_debugStatus = (cfiMemoryDebug); */


int read_paramaters(int parm_num, char escape, cfiScalar *params, FILE *fptr)
{
  int i=0,j,buf_pos=0;
  char buffer[MAX_STR_LEN],parm[MAX_STR_LEN];
  char *str_pos;

  while (i<parm_num)
  {
    if(fgets (buffer , MAX_STR_LEN , fptr)==NULL)
    	return i;

    str_pos = strchr(buffer,escape);
    /*printf ("string:  %s",buffer);*/
    if((str_pos-buffer+1 >= 0) && (str_pos-buffer+1 <= MAX_STR_LEN))
    {
	    strncpy (parm,buffer,str_pos-buffer+1);
	    parm[str_pos-buffer+1]='\0';
	    params[i] = atof(parm);
	    /*printf ("From: %s   we get: %f    ",parm,params[i]);*/
	    i++;
    }
    /*printf ("--> found at %d \n",str_pos-buffer+1);*/
  }
  return i;
}



main(int argc, char **argv)
{
	cfiErrorStatus errorStatus;

	FILE *file_ptr;
	char  fileName[MAX_STR_LEN];
	char ATN_FILE_NAME[MAX_STR_LEN];
	char OUT_FILE_NAME[MAX_STR_LEN];

	/* DEBUG: on */
	/* changed from float */
	cfiScalar *cellData;

	cfiScalar *params;
	cfiCounter *cellIndex;
	cfiCounter nonZeroCells;
	cfiCounter dataWritten;
	cfiCounter numberToBeWritten,index;


	cfiMesh   *atnMesh;
	cfiArray  *atnMeshDataArray;
	cfiScalar *atnMeshValue;

	cfiMesh   *rayMesh;
	cfiArray  *rayMeshDataArray;
	cfiVector *rayVector;
	cfiScalar *rayValue;

	cfiVector *atnVector;
	cfiScalar *atnVectorValue;

	cfiMatrix *F;
	cfiScalar *FValue;

	cfiMatrix *sampleIndices;
	cfiScalar *sampleIndexValue;

	cfiCounter indices[3];
	cfiScalar  origin[3];
	cfiScalar  coords[3];

	cfiScalar  voxel_min_width;
	cfiScalar  fov_radius;

	cfiScalar  detector_origin_x, detector_origin_y, detector_origin_z;
	cfiScalar  detector_dx, detector_dy, detector_dz;
	cfiScalar  detector_x, detector_y, detector_z;

	cfiScalar  ptresp_cone_radius;
	cfiScalar  ptresp_bin_width;
	cfiCounter ptresp_count;

	cfiScalar  **aperture_factor;
	cfiScalar  aperture_factor_sum, aperture_factor_sum2;
	cfiScalar  rT;

	cfiScalar  ptresp_origin_x, ptresp_origin_y, ptresp_origin_z;
	cfiScalar  ptresp_dx, ptresp_dy, ptresp_dz;
	cfiScalar  ptresp_x, ptresp_y, ptresp_z;

	cfiScalar  ray_dx, ray_dy, ray_dz, ray_ds, ray_sample_ds;
	cfiScalar  radial_dx, radial_dy, radial_dz;
	cfiScalar  tangential_dx, tangential_dy, tangential_dz;

	cfiCounter atn_voxel_count, voxel_count, voxel_xy_count;
	cfiCounter proj_count, ray_sample_count;
	cfiCounter angle, drow, dcol, prow, pcol;
	cfiScalar *atn_ptr, *ray_ptr, *end_ray;
	cfiScalar  atn_val;

	cfiScalar *Frow_ptr, *index_ptr;
	cfiCounter Fcol0, Fcol1, Fcol2, Fcol3;
	cfiCounter Fcol4, Fcol5, Fcol6, Fcol7;
	cfiScalar  u0, u1, u2, om_u0, om_u1, om_u2;
	cfiCounter atn_zpad, img_zmax;

	cfiCounter final_scafac;
	errorStatus = cfiSuccess;

	/* DEBUG: on */
	printf( "cfiCounter size: %ld bytes\n", sizeof( cfiCounter ) );
	printf( "cfiScalar size: %ld bytes\n", sizeof( cfiScalar ) );

	if((argc != 3) && (argc != 4))
	{
		fprintf(stderr, "Usage: > %s sysmat.fname params.fname (attn_matrix.fname)\n",argv[0]);
		return 1;
	}
	if( (file_ptr=fopen(argv[2], "rb"))==NULL)
	{
		fprintf(stderr, "Failed to open parameters file %s\n", argv[2]);
		exit(1);
	}

	params = (cfiScalar *) malloc(sizeof(cfiScalar)*PARM_NUM);
	if(index=read_paramaters(PARM_NUM, ':', params, file_ptr) != PARM_NUM)
	{
		fprintf(stderr, "Expected %d paramaters in file %s, only read %d paramaters\n",PARM_NUM, argv[2],index);
		exit(1);
	}
	fclose(file_ptr);
	strcpy(OUT_FILE_NAME,argv[1]);

	cfiCounter ANGLE_COUNT = (cfiCounter)(params[0]);
	cfiScalar ANGLE_START = params[1]/360;
	cfiScalar ANGLE_MODIFIER = params[2]/360;
	cfiScalar ANGLE_DIRECTION = params[3];

	cfiCounter XY_DIRECTION = (cfiCounter)(params[4]);
	cfiCounter Z_DIRECTION = (cfiCounter)(params[5]);
	cfiCounter DETECTOR_XY_COUNT = (cfiCounter)(params[6]);
	cfiCounter DETECTOR_Z_COUNT = (cfiCounter)(params[7]);

	cfiCounter VOXEL_X_COUNT = (cfiCounter)(params[8]);
	cfiCounter VOXEL_Y_COUNT = (cfiCounter)(params[9]);
	cfiCounter VOXEL_Z_COUNT = (cfiCounter)(params[10]);
	cfiCounter ATN_Z_COUNT = (cfiCounter)(params[11]);

	cfiScalar VOXEL_X_WIDTH = params[12];
	cfiScalar VOXEL_Y_WIDTH = params[13];
	cfiScalar VOXEL_Z_WIDTH = params[14];
	cfiScalar DETECTOR_BIN_WIDTH = params[15];

	cfiScalar DETECTOR_RADIUS = params[16];
	cfiScalar FOCAL_LENGTH = params[17];

	cfiScalar COLLIMATOR_RADIUS = params[18];
	cfiScalar COLLIMATOR_LENGTH = params[19];
	cfiScalar COLLIMATOR_OFFSET = params[20];

	cfiCounter RAY_SAMPLES_PER_VOXEL = (cfiCounter)(params[21]);
	cfiCounter RAYS_PER_VOXEL = (cfiCounter)(params[22]);
	cfiCounter MAX_PTRESP_COUNT = (cfiCounter)(params[23]);

	fprintf(stderr, "angle %d  %f  %f  %f\n",ANGLE_COUNT,ANGLE_START,ANGLE_DIRECTION,ANGLE_MODIFIER);
	fprintf(stderr, "detect  %d  %d  %d  %d\n",XY_DIRECTION,Z_DIRECTION,DETECTOR_XY_COUNT,DETECTOR_Z_COUNT);
	fprintf(stderr, "counts  %d  %d  %d  %d\n",VOXEL_X_COUNT,VOXEL_Y_COUNT,VOXEL_Z_COUNT,ATN_Z_COUNT);
	fprintf(stderr, "width  %f  %f  %f  %f\n",VOXEL_X_WIDTH,VOXEL_Y_WIDTH,VOXEL_Z_WIDTH,DETECTOR_BIN_WIDTH);
	fprintf(stderr, "rads  %f  %f\n",DETECTOR_RADIUS,FOCAL_LENGTH);
	fprintf(stderr, "colm  %f  %f  %f    rays  %d %d %d \n",COLLIMATOR_RADIUS,COLLIMATOR_LENGTH,COLLIMATOR_OFFSET,RAY_SAMPLES_PER_VOXEL,RAYS_PER_VOXEL,MAX_PTRESP_COUNT);	/*  allocate and read the attenuation map */

	aperture_factor = (cfiScalar **)   malloc(sizeof(cfiScalar *)  *  MAX_PTRESP_COUNT);
	for(index = 0; index<MAX_PTRESP_COUNT; index++)
		aperture_factor[index] = (cfiScalar *)   malloc(sizeof(cfiScalar)  *  MAX_PTRESP_COUNT);

	indices[xAxis] = VOXEL_X_COUNT;
	indices[yAxis] = VOXEL_Y_COUNT;
	indices[zAxis] = ATN_Z_COUNT;
	atn_voxel_count = VOXEL_X_COUNT * VOXEL_Y_COUNT * ATN_Z_COUNT;

	origin[xAxis] = (-((cfiScalar)(VOXEL_X_COUNT - 1)) / 2.0) * VOXEL_X_WIDTH;
	origin[yAxis] = (-((cfiScalar)(VOXEL_Y_COUNT - 1)) / 2.0) * VOXEL_Y_WIDTH;
	origin[zAxis] = (-((cfiScalar)(ATN_Z_COUNT   - 1)) / 2.0) * VOXEL_Z_WIDTH;

	/* DEBUG: on */
	printf( "Indices: ( %d, %d, %d )\n", indices[ xAxis ], indices[ yAxis ], indices[ zAxis ] );
	printf( "Origin: ( %f, %f, %f )\n", origin[ xAxis ], origin[ yAxis ], origin[ zAxis ] );

	if ( ( ( atnMesh = allocateMesh( ( cfiMeshType ) uniformQuadCornerType, 3, ( cfiDataType ) scalarType, 3, indices ) )
		== ( cfiMesh *) NULL )
	||   ( ( atnMeshDataArray = getMeshDataArrayPtr( atnMesh ) )
		== ( cfiArray * ) NULL )
	||   ( ( atnMeshValue = ( cfiScalar * )getArrayValuePtr( atnMeshDataArray, ( cfiCounter * )NULL ) )
		== (cfiScalar *)NULL ) )
	{
		FAILURE_MSG( "cannot allocate attenuation map ");
		( void )( exit( cfiFailure ) );
	}

	errorStatus |= setUniQuadMeshOrigin(origin, atnMesh);
	coords[xAxis] = VOXEL_X_WIDTH;
	coords[yAxis] = 0.0;
	coords[zAxis] = 0.0;
	errorStatus |= setUniQuadMeshAxis(coords, atnMesh, xAxis);
	coords[xAxis] = 0.0;
	coords[yAxis] = VOXEL_Y_WIDTH;
	coords[zAxis] = 0.0;
	errorStatus |= setUniQuadMeshAxis(coords, atnMesh, yAxis);
	coords[xAxis] = 0.0;
	coords[yAxis] = 0.0;
	coords[zAxis] = VOXEL_Z_WIDTH;
	errorStatus |= setUniQuadMeshAxis(coords, atnMesh, zAxis);


	if(argc == 4)
	{
		strcpy(ATN_FILE_NAME,argv[3]);
		fprintf(stdout, "reading %s\n", ATN_FILE_NAME);

		if ( ( (file_ptr = (FILE *)fopen(ATN_FILE_NAME, "r"))== (FILE *)NULL )
				||   ( fread(atnMeshValue, sizeof(cfiScalar), atn_voxel_count, file_ptr) != atn_voxel_count )
				||   ( fclose(file_ptr) == EOF ) )
		{
			FAILURE_MSG("cannot read file \n");
			(void)exit(cfiFailure);
		}
	}
	else
		/*memset(atnMeshValue, 0, atn_voxel_count * sizeof(cfiScalar));*/
		for(index = 0; index<atn_voxel_count; index++)
			atnMeshValue[index] = 1.0e-030;

	/* allocate projection ray samples and attenuation factors */

	voxel_min_width = VOXEL_X_WIDTH;
	if ( voxel_min_width > VOXEL_Y_WIDTH )
	{
		voxel_min_width = VOXEL_Y_WIDTH;
	}
	if ( voxel_min_width > VOXEL_Z_WIDTH )
	{
		voxel_min_width = VOXEL_Z_WIDTH;
	}

	fov_radius = DETECTOR_RADIUS - COLLIMATOR_OFFSET - COLLIMATOR_LENGTH;

	ptresp_cone_radius = 2.0 * COLLIMATOR_RADIUS * FOCAL_LENGTH / COLLIMATOR_LENGTH;
	ptresp_bin_width = (FOCAL_LENGTH / (DETECTOR_RADIUS + fov_radius)) * (voxel_min_width / RAYS_PER_VOXEL);
	ptresp_count = 1 + (2 * ((cfiCounter)(ptresp_cone_radius / ptresp_bin_width)));
	if ( ptresp_count > MAX_PTRESP_COUNT )
	{
		ptresp_count = MAX_PTRESP_COUNT;
	}

	ray_dx = DETECTOR_RADIUS + fov_radius;
	ray_dy = (((cfiScalar)(ptresp_count - 1)) / 2.0) * ptresp_bin_width;
	ray_dz = (((cfiScalar)(ptresp_count - 1)) / 2.0) * ptresp_bin_width;
	ray_ds = sqrt((ray_dx * ray_dx) + (ray_dy * ray_dy) + (ray_dz * ray_dz));
	ray_sample_ds = voxel_min_width / RAY_SAMPLES_PER_VOXEL;

	ray_sample_count = 1 + ((cfiCounter)(ray_ds / ray_sample_ds));
	indices[axis0] = ray_sample_count;

	fprintf(stderr, "ray_sample_count = %d\n", ray_sample_count);
	fprintf(stderr, "ptresp_count = %d\n", ptresp_count);

	if ( ( (rayMesh = allocateMesh((cfiMeshType)uniformQuadCornerType, 3, (cfiDataType)scalarType, 1, indices))
		== (cfiMesh *)NULL )
	||   ( (rayMeshDataArray = getMeshDataArrayPtr(rayMesh))
		== (cfiArray *)NULL )
	||   ( (rayValue = (cfiScalar *)getArrayValuePtr(rayMeshDataArray, (cfiCounter *)NULL))
		== (cfiScalar *)NULL )

	||   ( (rayVector = allocateVector((cfiDataType)scalarType, 0))
		== (cfiVector *)NULL )
	||   ( setVectorValueCount(indices[axis0], rayVector) != cfiSuccess )
	||   ( setVectorValuePtr(rayValue, rayVector) != cfiSuccess )

	||   ( (atnVector = allocateVector((cfiDataType)scalarType, indices[axis0]))
		== (cfiVector *)NULL )
	||   ( (atnVectorValue = (cfiScalar *)getVectorValuePtr(atnVector, 0))
		== (cfiScalar *)NULL ) )
	{
		FAILURE_MSG("cannot allocate rays");
		(void)exit(cfiFailure);
	}

	/*
	 * allocate an F matrix to map image space to projections
	 * for one angle, and allocate storage for sample indices
	 * from which F matrix can be calculated
	 */

	proj_count = DETECTOR_XY_COUNT * DETECTOR_Z_COUNT;
	voxel_xy_count = VOXEL_X_COUNT * VOXEL_Y_COUNT;
	voxel_count = VOXEL_X_COUNT * VOXEL_Y_COUNT * VOXEL_Z_COUNT;

	atn_zpad = (ATN_Z_COUNT - VOXEL_Z_COUNT) / 2;
	img_zmax = VOXEL_Z_COUNT - 1;

	if ( ( (F = allocateMatrix((cfiDataType)scalarType, 1, voxel_count))
		== (cfiMatrix *)NULL )
	||   ( (FValue = (cfiScalar *)getMatrixValuePtr(F, 0, 0))
		== (cfiScalar *)NULL )
	||   ( (sampleIndices = allocateMatrix((cfiDataType)scalarType, ray_sample_count, 3))
		== (cfiMatrix *)NULL )
	||   ( (sampleIndexValue = (cfiScalar *)getMatrixValuePtr(sampleIndices, 0, 0))
		== (cfiScalar *)NULL ) )
	{
		FAILURE_MSG("cannot allocate F matrix");
		(void)exit(cfiFailure);
	}

	/*
	 * calculate the aperture autocorrelation factors
	 */

	aperture_factor_sum = 0.0;
	for ( prow = 0; prow < ptresp_count; prow++ )
	{
	    for ( pcol = 0; pcol < ptresp_count; pcol++ )
	    {
		rT = sqrt(((prow - (ptresp_count / 2)) * (prow - (ptresp_count / 2))) + ((pcol - (ptresp_count / 2)) * (pcol - (ptresp_count / 2)))) * ptresp_bin_width;
		rT *= COLLIMATOR_LENGTH / FOCAL_LENGTH;
		rT /= (2.0 * COLLIMATOR_RADIUS);
		if ( rT >= 1.0 )
		{
			aperture_factor[prow][pcol] = 0.0;
		}
		else
		{
			aperture_factor[prow][pcol] = acos(rT) - (rT * sqrt(1.0 - (rT * rT)));
		}
		aperture_factor_sum += aperture_factor[prow][pcol];
	    }
	}
	aperture_factor_sum2 = 0.0;
	for ( prow = 0; prow < ptresp_count; prow++ )
	{
	    for ( pcol = 0; pcol < ptresp_count; pcol++ )
	    {
		aperture_factor[prow][pcol] /= aperture_factor_sum;
		aperture_factor_sum2 += aperture_factor[prow][pcol];
	    }
	}
	fprintf(stderr, "aperture factors (sum = %f)\n", aperture_factor_sum2);
	for ( prow = 0; prow <= ptresp_count / 2; prow++ )
	{
	    for ( pcol = 0; pcol <= ptresp_count / 2; pcol++ )
	    {
		fprintf(stderr, "%2d %2d %f\n", prow, pcol, aperture_factor[prow][pcol]);
	    }
	}

	/*
	 * calculate the F matrix
	 */
#ifdef SPARSE_FORMAT

	  cellData = ( float * ) malloc( ( size_t ) sizeof( float ) *  VOXEL_X_COUNT * VOXEL_Y_COUNT * VOXEL_Z_COUNT );
	  cellIndex = ( cfiCounter * ) malloc( ( size_t ) sizeof( cfiCounter ) *  VOXEL_X_COUNT * VOXEL_Y_COUNT * VOXEL_Z_COUNT );

	  fprintf( stderr, "writing %s \n", OUT_FILE_NAME );

	  if ( ( file_ptr = (FILE *)fopen( OUT_FILE_NAME, "w" ) ) == ( FILE * )NULL )
	  {
		FAILURE_MSG("cannot write file");
		(void)exit(cfiFailure);
	  }

        numberToBeWritten=ANGLE_COUNT*DETECTOR_XY_COUNT*DETECTOR_Z_COUNT;
        /* DEBUG: on */
        printf( "Writing %d in %d bytes...", numberToBeWritten, ( int ) sizeof(cfiCounter) );
        fwrite(&numberToBeWritten , sizeof(cfiCounter), 1, file_ptr);


        numberToBeWritten=VOXEL_X_COUNT*VOXEL_Y_COUNT*VOXEL_Z_COUNT;
        /* DEBUG: on */
        printf( "Writing %d in %d bytes...", numberToBeWritten, ( int ) sizeof(cfiCounter) );
        fwrite(&numberToBeWritten , sizeof(cfiCounter), 1, file_ptr);
#endif

	/*for ( angle = 0; angle < ANGLE_COUNT; angle++ )*/
	/*for ( angle = ANGLE_COUNT-1; angle >=0; angle-- )*/
	for ( angle = 0; angle < ANGLE_COUNT; angle++ )
	{
#ifndef SPARSE_FORMAT
	  (void)sprintf(fileName, "%s.%03d",OUT_FILE_NAME, angle);
	  fprintf(stderr, "writing %s ", fileName);

	  if ( (file_ptr = (FILE *)fopen(fileName, "w")) == (FILE *)NULL )
	  {
		FAILURE_MSG("cannot write file");
		(void)exit(cfiFailure);
	  }
#endif
	  /*
	   * set up geometry of detector plane and
	   * plane containing base of point response cone
	   */

	  radial_dx = cos(2.0 * M_PI * (ANGLE_START + ((cfiScalar)(angle*ANGLE_MODIFIER*ANGLE_DIRECTION)) / ((cfiScalar)ANGLE_COUNT)));
	  radial_dy = sin(2.0 * M_PI * (ANGLE_START + ((cfiScalar)(angle*ANGLE_MODIFIER*ANGLE_DIRECTION)) / ((cfiScalar)ANGLE_COUNT)));
	  radial_dz = 0.0;

	  fprintf(stderr,"dx: %.2f  dy: %.2f .an: %d | ",radial_dx,radial_dy,angle);

	  tangential_dx = (-radial_dy);
	  tangential_dy =   radial_dx;
	  tangential_dz =   0.0;

	  detector_dx = DETECTOR_BIN_WIDTH * tangential_dx;
	  detector_dy = DETECTOR_BIN_WIDTH * tangential_dy;
	  detector_dz = DETECTOR_BIN_WIDTH;

	  ptresp_dx = ptresp_bin_width * (-tangential_dx);
	  ptresp_dy = ptresp_bin_width * (-tangential_dy);
	  ptresp_dz = ptresp_bin_width;

	  detector_origin_x = (DETECTOR_RADIUS * radial_dx) + (((cfiScalar)(DETECTOR_XY_COUNT - 1) / 2.0) * (-detector_dx));
	  detector_origin_y = (DETECTOR_RADIUS * radial_dy) + (((cfiScalar)(DETECTOR_XY_COUNT - 1) / 2.0) * (-detector_dy));
	  detector_origin_z = (((cfiScalar)(DETECTOR_Z_COUNT - 1)) / 2.0) * (-detector_dz);

	  fprintf(stderr,"x: %.1f  y: %.1f z: %.1f | ",detector_origin_x,detector_origin_y,detector_origin_z);

	  for ( drow = 0; drow < DETECTOR_Z_COUNT; drow++ )
	  {
	    fprintf(stderr, ".");

	    if(Z_DIRECTION)
		    detector_z = detector_origin_z + (drow * detector_dz);
	    else
		    detector_z = detector_origin_z + ((DETECTOR_Z_COUNT-drow-1) * detector_dz);


	    origin[zAxis] = detector_z;

	    ptresp_origin_z = detector_z + (((cfiScalar)(ptresp_count - 1)) / 2.0) * (-ptresp_dz);

	    for ( dcol = 0; dcol < DETECTOR_XY_COUNT; dcol++ )
	    {
		errorStatus |= zeroMatrix(F);
		Frow_ptr = FValue;
	      /*
	       * calculate position on detector plane
	       */

	       if(XY_DIRECTION)
	       {
		       detector_x = detector_origin_x + (dcol * detector_dx);
		       detector_y = detector_origin_y + (dcol * detector_dy);
	       }
	       else
	       {
		       detector_x = detector_origin_x + ((DETECTOR_XY_COUNT-dcol-1) * detector_dx);
		       detector_y = detector_origin_y + ((DETECTOR_XY_COUNT-dcol-1) * detector_dy);
	       }
	      origin[xAxis] = detector_x;
	      origin[yAxis] = detector_y;
	      errorStatus |= setUniQuadMeshOrigin(origin, rayMesh);

	      ptresp_origin_x = detector_x + (FOCAL_LENGTH * (-radial_dx)) + (((cfiScalar)(ptresp_count - 1) / 2.0) * (-ptresp_dx));
	      ptresp_origin_y = detector_y + (FOCAL_LENGTH * (-radial_dy)) + (((cfiScalar)(ptresp_count - 1) / 2.0) * (-ptresp_dy));

	      for ( prow = 0; prow < ptresp_count; prow++ )
	      {
		ptresp_z = ptresp_origin_z + (prow * ptresp_dz);
		ray_dz = ptresp_z - detector_z;

		for ( pcol = 0; pcol < ptresp_count; pcol++ )
		{
		  /*
		   * calculate position on plane containing base
		   * of point response cone and set up ray geometry
		   */

		  ptresp_x = ptresp_origin_x + (pcol * ptresp_dx);
		  ptresp_y = ptresp_origin_y + (pcol * ptresp_dy);

		  ray_dx = ptresp_x - detector_x;
		  ray_dy = ptresp_y - detector_y;
		  ray_ds = sqrt((ray_dx * ray_dx) + (ray_dy * ray_dy) + (ray_dz * ray_dz)) / ray_sample_ds;

		  coords[xAxis] = ray_dx / ray_ds;
		  coords[yAxis] = ray_dy / ray_ds;
		  coords[zAxis] = ray_dz / ray_ds;
		  errorStatus |= setUniQuadMeshAxis(coords, rayMesh, (cfiDataAxis)axis0);

		  /*
		   * calculate attenuation along the ray,
		   * and gather information in sampleIndices array
		   * so that F matrix can be calculated
		   *
		   * treat attenuation as being piecewise constant
		   * in segments centered at sampling points
		   *
		   * treat the activity in each segment as being concentrated
		   * at the center of the segment (i.e., an impulse)
		   *
		   * thus, the activity in each segment is attenuated
		   * over one-half of the segment in which it lies, as well
		   * over subsequent segments leading to the detector
		   *
		   * NOTE:  the vector rayVector and the pointer rayValue
		   * provide access to the values in rayMesh
		   */

		  errorStatus |= mapUniQuadMeshValues(atnMesh, (cfiMeshSamplingMethod)linearInterp, rayMesh, sampleIndices);
		  errorStatus |= scaleVector(ray_sample_ds, rayVector);

		  atn_ptr = atnVectorValue;
		  ray_ptr = rayValue;
		  end_ray = ray_ptr + ray_sample_count - 1;

		  (*atn_ptr) = 0.0;
		  for ( ; ray_ptr < end_ray; ray_ptr++, atn_ptr++ )
		  {
		    (*(atn_ptr + 1)) = (*atn_ptr) + (*ray_ptr);
		    (*atn_ptr) += 0.5 * (*ray_ptr);
		    (*atn_ptr) = exp(-(*atn_ptr));
		  }
		  (*atn_ptr) += 0.5 * (*ray_ptr);
		  (*atn_ptr) = exp(-(*atn_ptr));

		  /*
		   * update row in F matrix corresponding
		   * to this detector element at this angle
		   */

		  index_ptr = sampleIndexValue;
		  atn_ptr = atnVectorValue;
		  ray_ptr = rayValue;
		  end_ray = ray_ptr + ray_sample_count;

		  for ( ; ray_ptr < end_ray; ray_ptr++, atn_ptr++ )
		  {
		     if ( (*ray_ptr) > 0.0 )
		     {
			/*
			 * adjust for z-axis padding in attenuation map
			 */
			index_ptr[2] -= atn_zpad;

			if ( index_ptr[2] < 0.0 )
			{
				/*
				 * make sure upper 4 voxels calculated next
				 * are on bottom edge of image volume
				 */
				index_ptr[2] = -0.5;
			}
			else if (index_ptr[2] > (cfiScalar)img_zmax)
			{
				/*
				 * make sure lower 4 voxels calculated next
				 * are on top edge of image volume
				 */
				index_ptr[2] = (cfiScalar)img_zmax + 0.5;
			}

			/*
			 * calculate pointer offsets to F matrix elements
			 * corresponding to 8 voxels surrounding sample
			 */
			Fcol0 = ((cfiCounter)index_ptr[0]) + (((cfiCounter)index_ptr[1]) * VOXEL_X_COUNT) + (((cfiCounter)floor(index_ptr[2])) * voxel_xy_count);

			Fcol1 = Fcol0 + 1;

			Fcol2 = Fcol0 + VOXEL_X_COUNT;
			Fcol3 = Fcol1 + VOXEL_X_COUNT;

			Fcol4 = Fcol0 + voxel_xy_count;
			Fcol5 = Fcol1 + voxel_xy_count;
			Fcol6 = Fcol2 + voxel_xy_count;
			Fcol7 = Fcol3 + voxel_xy_count;

			u0 = index_ptr[0] - ((cfiCounter)index_ptr[0]);
			u1 = index_ptr[1] - ((cfiCounter)index_ptr[1]);
			u2 = index_ptr[2] - ((cfiCounter)floor(index_ptr[2]));

			om_u0 = 1.0 - u0;
			om_u1 = 1.0 - u1;
			om_u2 = 1.0 - u2;

			atn_val = (*atn_ptr) * aperture_factor[prow][pcol];

			if ( index_ptr[2] >= 0.0 )
			{
			   *(Frow_ptr + Fcol0) += atn_val * om_u0 * om_u1 * om_u2;
			   *(Frow_ptr + Fcol1) += atn_val *    u0 * om_u1 * om_u2;
			   *(Frow_ptr + Fcol2) += atn_val * om_u0 *    u1 * om_u2;
			   *(Frow_ptr + Fcol3) += atn_val *    u0 *    u1 * om_u2;
			}
			else
			{
			   /*
			    * in the attenuator, but below the image volume;
			    * thus, assign all of the contribution to four
			    * nearest voxels on bottom edge of image volume
			    * by executing these statements and the statements
			    * in top half of next if-else statement
			    */
			   *(Frow_ptr + Fcol4) += atn_val * om_u0 * om_u1 * om_u2;
			   *(Frow_ptr + Fcol5) += atn_val *    u0 * om_u1 * om_u2;
			   *(Frow_ptr + Fcol6) += atn_val * om_u0 *    u1 * om_u2;
			   *(Frow_ptr + Fcol7) += atn_val *    u0 *    u1 * om_u2;
			}

			if ( index_ptr[2] < (cfiScalar)img_zmax )
			{
			   *(Frow_ptr + Fcol4) += atn_val * om_u0 * om_u1 *    u2;
			   *(Frow_ptr + Fcol5) += atn_val *    u0 * om_u1 *    u2;
			   *(Frow_ptr + Fcol6) += atn_val * om_u0 *    u1 *    u2;
			   *(Frow_ptr + Fcol7) += atn_val *    u0 *    u1 *    u2;
			}
			else
			{
			   /*
			    * in the attenuator, but above the image volume;
			    * thus, assign all of the contribution to four
			    * nearest voxels on top edge of image volume
			    * by executing these statements and the statements
			    * in top half of previous if-else statement
			    */
			   *(Frow_ptr + Fcol0) += atn_val * om_u0 * om_u1 *    u2;
			   *(Frow_ptr + Fcol1) += atn_val *    u0 * om_u1 *    u2;
			   *(Frow_ptr + Fcol2) += atn_val * om_u0 *    u1 *    u2;
			   *(Frow_ptr + Fcol3) += atn_val *    u0 *    u1 *    u2;
			}

		     } /* end of if ( (*ray_ptr) > 0.0 ) */

		     index_ptr += 3;

		  } /* end of for ( ; ray_ptr < end_ray; ... ) loop */

		} /* end of for ( pcol = 0; ... ) loop */

	      } /* end of for ( prow = 0; ... ) loop */

	      /*Frow_ptr += voxel_count;*/

	  errorStatus |= scaleMatrix(1.0 / RAY_SAMPLES_PER_VOXEL, F);

#ifndef SPARSE_FORMAT
	  if ( (errorStatus & cfiFailure) )
	  {
		FAILURE_MSG("cannot project");
		(void)exit(cfiFailure);
	  }

	  if ( fwrite(FValue, sizeof(cfiScalar), voxel_count, file_ptr) != ( voxel_count) )
	  {
		FAILURE_MSG("cannot write file");
		(void)exit(cfiFailure);
	  }
#else
	      nonZeroCells = 0;
	      dataWritten = 0;

	      for ( index = 0; index < voxel_count ; index++ )
		{
		  if(FValue[index]!=0)
		    {
		      cellData[nonZeroCells] = (float)FValue[index];
		      cellIndex[nonZeroCells] = index;
	 	      nonZeroCells++;
		    }
		}


	      numberToBeWritten = dcol + drow*DETECTOR_XY_COUNT + angle * DETECTOR_Z_COUNT*DETECTOR_XY_COUNT;
	      dataWritten += fwrite(&numberToBeWritten, sizeof(cfiCounter), 1, file_ptr);
	      dataWritten += fwrite(&nonZeroCells, sizeof(cfiCounter), 1, file_ptr);
	      if(nonZeroCells)
		{
		  dataWritten += fwrite(cellIndex,sizeof(cfiCounter), nonZeroCells, file_ptr);
		  dataWritten += fwrite(cellData, sizeof(float), nonZeroCells, file_ptr);
		}

	      if ( dataWritten != (2 + 2 * nonZeroCells) )
		{
		  FAILURE_MSG("cannot write file.");
		  (void)exit(cfiFailure);
		}
#endif

	    } /* end of for ( dcol = 0; ... ) loop */

	  } /* end of for ( drow = 0; ... ) loop */

#ifndef SPARSE_FORMAT
	   if ( fclose(file_ptr) == EOF )
	  {
		FAILURE_MSG("cannot write file");
		(void)exit(cfiFailure);
	  }
#endif

	  fprintf(stderr, "\n");

	} /* end of for ( angle = 0; ... ) loop */

#ifdef SPARSE_FORMAT
	fclose(file_ptr);
	free(cellData);
	free(cellIndex);
#endif

	for(index = 0; index<MAX_PTRESP_COUNT; index++)
		free(aperture_factor[index]);
	free(aperture_factor);

} /* end of main() */
