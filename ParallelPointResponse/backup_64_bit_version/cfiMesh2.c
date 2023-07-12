/*
 * SCCS:  %Z%%M%  %I%  %G%  %U%
 */

/************************************************************************
 *									*
 * -----------------							*
 * FUNCTIONS HEREIN:							*
 * -----------------							*
 *									*
 * mapUniQuadMeshValues							*
 * --------------------							*
 *	Maps data values from one uniform quadrilateral mesh		*
 *	to another uniform quadrilateral mesh.  Linearly		*
 *	interpolates between the values at adjacent source		*
 *	mesh vertices, and stores the interpolated values		*
 *	in the destination mesh.					*
 *									*
 * getUniQuadMeshIndexTransform						*
 * ----------------------------						*
 *	Calculates the matrix that transforms floating-point		*
 *	indices	of one uniform quadrilateral mesh data array into	*
 *	floating-point indices of a second uniform quadrilateral	*
 *	mesh data array, such that the latter set of array indices	*
 *	corresponds to the same coordinate location as the former.	*
 *									*
 * sampleMeshValue							*
 * ---------------							*
 *	Samples a mesh data array at an arbitrary index location.	*
 *	Linearly interpolates between the values at the mesh vertices	*
 *	adjacent to the index location.					*
 *									*
 * allocateMesh								*
 * ------------								*
 *	Returns a pointer to a mesh of the specified geometry type,	*
 *	coordinate dimension, and data value type and dimensions.	*
 *	The mesh coordinates and data values are set to zero.		*
 *									*
 * deallocateMesh							*
 * --------------							*
 *	Deallocates memory for a mesh.					*
 *									*
 * writeMesh								*
 * ---------								*
 *	Writes a mesh to a stream.					*
 *									*
 * checkMesh								*
 * ---------								*
 *	Checks if a mesh has data values at				*
 *	and coordinates for its vertices.				*
 *									*
 * getUniQuadMeshOrigin							*
 * --------------------							*
 *	Returns the position vector for the origin of a uniform		*
 *	quadrilateral mesh, in a contiguous block of scalars.		*
 *									*
 * setUniQuadMeshOrigin							*
 * --------------------							*
 *	Sets the position vector for the origin of a uniform		*
 *	quadrilateral mesh, from a contiguous block of scalars.		*
 *									*
 * getUniQuadMeshAxis							*
 * ------------------							*
 *	Returns the direction vector for an axis of a uniform		*
 *	quadrilateral mesh, in a contiguous block of scalars.		*
 *									*
 * setUniQuadMeshAxis							*
 * ------------------							*
 *	Sets the direction vector for an axis of a uniform		*
 *	quadrilateral mesh, from a contiguous block of scalars.		*
 *									*
 * getMeshType								*
 * -----------								*
 *	Returns the type of geometry for a mesh.			*
 *									*
 * getMeshCoordDimension						*
 * ---------------------						*
 *	Returns the vertex coordinate dimension for a mesh.		*
 *									*
 * getMeshAxisCount							*
 * ----------------							*
 *	Returns the number of axes spanning a mesh.			*
 *									*
 * getMeshDataArrayPtr							*
 * -------------------							*
 *	Returns a pointer to the data array for a mesh.			*
 *									*
 ************************************************************************/



#include <stdio.h>
#include <math.h>

/*
 * "old"
 *
#include <malloc.h>
 *
 */

/*
 * for MacOS X
 */
#include <stdlib.h>

#include "cfiError.h"
#include "cfiTypes.h"
#include "cfiVector.h"
#include "cfiMatrix.h"
#include "cfiArray.h"
#include "cfiMesh2.h"



static char l_msg[MSG_LENGTH];		/* for compound messages */



/************************************************************************
 *									*
 * mapUniQuadMeshValues							*
 * --------------------							*
 *	Maps data values from one uniform quadrilateral mesh		*
 *	to another uniform quadrilateral mesh.  Linearly		*
 *	interpolates between the values at adjacent source		*
 *	mesh vertices, and stores the interpolated values		*
 *	in the destination mesh.					*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
mapUniQuadMeshValues
(
	const cfiMesh               *const fromMesh,	/* source mesh       */
	const cfiMeshSamplingMethod  samplingMethod,	/* sampling method   */
	cfiMesh                     *const toMesh,	/* destination mesh  */
	cfiMatrix                   *const sampleIndices/* sample indices    */
)

{
	cfiErrorStatus errorStatus;	/* the return value                  */

	cfiMeshStatus meshStatus;	/* status of the meshes              */

	cfiMatrixStatus  sampleIndexStatus;	/* status of index storage   */
	cfiScalar       *sampleIndex_ptr;

	cfiCounter fromAxisCount;	/* # axes spanning source mesh       */
	cfiCounter toAxisCount;		/* # axes spanning destination mesh  */

	const cfiArray *fromDataArray;	/* source mesh data array            */
	cfiDataType     fromValueType;	/* type of array values              */
	void           *fromValue_ptr;	/* pointer to array values           */

	cfiArray  *toDataArray;		/* destination mesh data array       */
	cfiMatrix *toValueMatrix;

	cfiCounter fromAxisDimensions[MAX_AXIS];	/* array dimensions  */
	cfiCounter toAxisDimensions[MAX_AXIS];

	cfiMatrix *indexTransform;	/* mesh data array index transform   */

	cfiScalar fromIndexIncrements[MAX_AXIS];	/* 1st row of above  */

	cfiVector *fromIndexVector;	/* mesh data array index storage     */
	cfiScalar *fromIndex_ptr;
	cfiVector *toIndexVector;
	cfiScalar *toIndex_ptr;

	cfiScalar *index_ptr;		/* source data array index pointers  */
	cfiScalar *indexIncr_ptr;
	cfiScalar *end_index;

	cfiCounter  neighborCount;	/* # of corners in source mesh voxel */
	cfiCounter *neighborOffsets;	/* array pointer offsets to corners  */
	cfiScalar  *neighborValues;	/* array values at corners           */

	cfiScalar *sampleValues;	/* buffers for sampled values        */
	cfiScalar *sample_ptr;
	cfiScalar *prevValues;
	cfiScalar *prev_ptr;
	cfiScalar *end_sample;

	cfiCounter fromAxis;		/* loop control                      */
	cfiCounter fromAxisOffset;
	cfiCounter lowerNeighbor;
	cfiCounter upperNeighbor;
	cfiCounter toAxis;
	cfiCounter toRow;
	cfiCounter toRowCount;
	cfiCounter toColCount;

	DEBUG_MSG(cfiTraceDebug, "mapping uniform quad-mesh values");

	meshStatus = (checkMesh(fromMesh) & checkMesh(toMesh));

	/*
	 * (sampleIndices matrix is checked later)
	 */

	if ( ( ! (meshStatus & meshHasValues) )
	||   ( ! (meshStatus & meshHasDimensions) )
	||   ( ! (meshStatus & meshHasCoords) )
	||   (    ( getMeshType(fromMesh) != uniformQuadCornerType )
	       && ( getMeshType(fromMesh) != uniformQuadCenterType ) )
	||   (    ( getMeshType(toMesh) != uniformQuadCornerType )
	       && ( getMeshType(toMesh) != uniformQuadCenterType ) )
	||   ( getMeshCoordDimension(fromMesh) != getMeshCoordDimension(toMesh) )
	||   ( (toAxisCount = getMeshAxisCount(toMesh))
		> (fromAxisCount = getMeshAxisCount(fromMesh)) )
	||   (    ( samplingMethod != (cfiMeshSamplingMethod)subNearestNeighbor )
	       && ( samplingMethod != (cfiMeshSamplingMethod)nearestNeighbor )
	       && ( samplingMethod != (cfiMeshSamplingMethod)linearInterp ) ) )
	{
			FAILURE_MSG("cannot map uniform quad-mesh values");
			return (cfiFailure);
	}

	/*
	 * allocate workspace
	 */

	errorStatus = cfiSuccess;

	toIndexVector = (cfiVector *)NULL;
	indexTransform = (cfiMatrix *)NULL;
	fromIndexVector = (cfiVector *)NULL;
	neighborOffsets = (cfiCounter *)NULL;
	neighborValues = (cfiScalar *)NULL;
	sampleValues = (cfiScalar *)NULL;
	prevValues = (cfiScalar *)NULL;
	toValueMatrix = (cfiMatrix *)NULL;

	neighborCount = (1L << fromAxisCount);

	if ( ( (toIndexVector = allocateVector((cfiDataType)scalarType, toAxisCount + 1))
		== (cfiVector *)NULL )
	||   ( (toIndex_ptr = (cfiScalar *)getVectorValuePtr(toIndexVector, 0))
		== (cfiScalar *)NULL )
	||   ( (indexTransform = allocateMatrix((cfiDataType)scalarType, toAxisCount + 1, fromAxisCount))
		== (cfiMatrix *)NULL )
	||   ( (fromIndexVector = allocateVector((cfiDataType)scalarType, fromAxisCount))
		== (cfiVector *)NULL )
	||   ( (fromIndex_ptr = (cfiScalar *)getVectorValuePtr(fromIndexVector, 0))
		== (cfiScalar *)NULL )
/* LINTED */
	||   ( (neighborOffsets = (cfiCounter *)malloc((size_t)(neighborCount * sizeof(cfiCounter))))
		== (cfiCounter *)NULL )
/* LINTED */
	||   ( (neighborValues = (cfiScalar *)malloc((size_t)(neighborCount * sizeof(cfiScalar))))
		== (cfiScalar *)NULL ) )
	{
		errorStatus |= cfiFailure;
		goto cleanupAndReturn;
	}

	end_index = fromIndex_ptr + fromAxisCount;

	/*
	 * get the transform that converts destination mesh data array
	 * indices to source mesh data array indices, and initialize
 	 * the destination mesh data array homogeneous index (see the
	 * comments in getUniQuadMeshIndexTransform() below)
	 */

	errorStatus |= getUniQuadMeshIndexTransform(toMesh, fromMesh, indexTransform);
/* LINTED */
	errorStatus |= getMatrixValues(indexTransform, 0, UNDEFINED_COUNT, fromIndexIncrements);

	toIndex_ptr[toAxisCount] = 1.0;		/* (other indices already 0) */

	if ( (errorStatus & cfiFailure) )
	{
		goto cleanupAndReturn;
	}

	/*
	 * set up faster access to mesh data arrays
	 */

	fromDataArray = getMeshDataArrayPtr(fromMesh);
	fromValueType = getArrayValueType(fromDataArray);
	fromValue_ptr = getArrayValuePtr(fromDataArray, (cfiCounter *)NULL);

	toDataArray = getMeshDataArrayPtr(toMesh);
	toAxisDimensions[0] = getArrayDimension(toDataArray, (cfiDataAxis)0);
	toColCount = toAxisDimensions[0];
	toRowCount = 1;
	for ( toAxis = 1; toAxis < toAxisCount; toAxis++ )
	{
		toAxisDimensions[toAxis] = getArrayDimension(toDataArray, (cfiDataAxis)toAxis);
		toRowCount *= toAxisDimensions[toAxis];
	}

/* LINTED */
	if ( ( (sampleValues = (cfiScalar *)malloc((size_t)(toColCount * sizeof(cfiScalar))))
		== (cfiScalar *)NULL )
/* LINTED */
	||   ( (prevValues = (cfiScalar *)malloc((size_t)(toColCount * sizeof(cfiScalar))))
		== (cfiScalar *)NULL )
	||   ( (toValueMatrix = allocateMatrix(getArrayValueType(toDataArray), 0, 0))
		== (cfiMatrix *)NULL ) )
	{
		errorStatus |= cfiFailure;
		goto cleanupAndReturn;
	}

	errorStatus |= setMatrixDimension((cfiDataAxis)rowAxis, toRowCount, toValueMatrix);
	errorStatus |= setMatrixDimension((cfiDataAxis)columnAxis, toColCount, toValueMatrix);
	errorStatus |= setMatrixValuePtr(getArrayValuePtr(toDataArray, (cfiCounter *)NULL), toValueMatrix);

	/*
	 * calculate the pointer offsets between neighboring columns,
	 * rows, planes, etc., in the source mesh data array (see the
	 * comments in sampleMeshValue() below for a diagram showing
	 * the layout and indexing of neighboring values in 3-D)
	 */

	fromAxisOffset = 1;		/* pointer offset between columns */

	neighborOffsets[0] = 0;
	neighborOffsets[1] = 1;

	fromAxisDimensions[0] = getArrayDimension(fromDataArray, (cfiDataAxis)0);

	for ( fromAxis = 1; fromAxis < fromAxisCount; fromAxis++ )
	{
		fromAxisOffset *= fromAxisDimensions[fromAxis - 1];

		upperNeighbor = (1L << fromAxis);
		for ( lowerNeighbor = 0; lowerNeighbor < upperNeighbor; lowerNeighbor++ )
		{
			neighborOffsets[upperNeighbor + lowerNeighbor] = neighborOffsets[lowerNeighbor] + fromAxisOffset;
		}

		fromAxisDimensions[fromAxis] = getArrayDimension(fromDataArray, (cfiDataAxis)fromAxis);
	}

	/*
	 * process the destination mesh vertices row by row,
	 * after checking optional storage for sample indices
	 */

	sampleIndexStatus = checkMatrix(sampleIndices);

	if ( (sampleIndexStatus & matrixHasValues) )
	{
		if ( ( ! (sampleIndexStatus & matrixHasDimensions) )
		||   ( getMatrixDimension(sampleIndices, (cfiDataAxis)columnAxis) != fromAxisCount )
		||   ( getMatrixDimension(sampleIndices, (cfiDataAxis)rowAxis) != (toRowCount * toColCount) )
		||   ( getMatrixValueType(sampleIndices) != scalarType )
		||   ( (sampleIndex_ptr = (cfiScalar *)getMatrixValuePtr(sampleIndices, 0, 0))
			== (cfiScalar *)NULL ) )
		{
			FAILURE_MSG("cannot store sample indices");
			return (cfiFailure);
		}
	}
	else
	{
		sampleIndex_ptr = (cfiScalar *)NULL;
	}

	toRow = 0;
	while ( toRow < toRowCount )
	{
		/*
		 * transform from the data array indices for the first
		 * vertex in the destination mesh row, to source mesh
		 * data array indices
		 */

		errorStatus |= vectorTimesMatrix(toIndexVector, indexTransform, fromIndexVector);

		if ( ( samplingMethod == (cfiMeshSamplingMethod)subNearestNeighbor )
		||   ( samplingMethod == (cfiMeshSamplingMethod)nearestNeighbor ) )
		{
			/*
			 * prepare for rounding the scalar indices
			 * to the nearest integers via the index
			 * truncations in sampleMeshValue()
			 */

			index_ptr = fromIndex_ptr;
			while ( index_ptr < end_index )
			{
				/*
				 * note that the function sampleMeshValue()
				 * below returns zero for any non-positive
				 * index--therefore add a half in all cases,
				 * so that positive indices are rounded
				 * correctly (if a half were subtracted to
				 * round an index that is negative initially
				 * and is positive later as one increments
				 * along the row, then the rounding would be
				 * incorrect for the positive values)
				 */

				(*index_ptr++) += 0.5;
			}
		}

		if ( (errorStatus & cfiFailure) )
		{
			goto cleanupAndReturn;
		}

		/*
		 * process a row of destination mesh vertices
		 */

		sample_ptr = sampleValues;

		end_sample = sampleValues + toColCount;
		while ( sample_ptr < end_sample )
		{
			/*
			 * sample the source data array at
			 * the destination mesh vertex location
			 */

			sampleMeshValue(samplingMethod, fromAxisCount, fromAxisDimensions, fromIndex_ptr, fromValueType, fromValue_ptr, neighborCount, neighborOffsets, neighborValues, sample_ptr);

			if ( samplingMethod == (cfiMeshSamplingMethod)subNearestNeighbor )
			{
				/*
				 * interpolate the previously skipped sample
				 * (sample skipped at end of row will be 0)
				 */

				if ( sample_ptr > sampleValues )
				{
					*(sample_ptr - 1) = ((*sample_ptr) + (*(sample_ptr - 2))) / 2.0;
				}

				/*
				 * skip the next sample
				 */

				indexIncr_ptr = fromIndexIncrements;
				index_ptr = fromIndex_ptr;
				while ( index_ptr < end_index )
				{
					(*index_ptr++) += (*indexIncr_ptr++);
				}

				sample_ptr++;
			}
			else if ( sampleIndex_ptr != (cfiScalar *)NULL )
			{
				/*
				 * store sample indices
				 */

				index_ptr = fromIndex_ptr;
				while ( index_ptr < end_index )
				{
					(*sampleIndex_ptr++) = (*index_ptr++);
				}
			}

			/*
			 * increment the source data array indices, to move
			 * to the next destination mesh vertex location
			 */

			indexIncr_ptr = fromIndexIncrements;
			index_ptr = fromIndex_ptr;
			while ( index_ptr < end_index )
			{
				(*index_ptr++) += (*indexIncr_ptr++);
			}

			sample_ptr++;

		} /* end of while ( sample_ptr < end_sample ) loop */

		/*
		 * store the row of samples in the destination data array
		 */

/* LINTED */
		errorStatus |= setMatrixValues(sampleValues, toValueMatrix, toRow, UNDEFINED_COUNT);

		if ( samplingMethod == (cfiMeshSamplingMethod)subNearestNeighbor )
		{
			/*
			 * interpolate the previously skipped row
			 */

			if ( toIndex_ptr[1] > 0.0 )
			{
				sample_ptr = sampleValues;
				prev_ptr = prevValues;

				end_sample = sampleValues + toColCount;
				while ( sample_ptr < end_sample )
				{
					(*prev_ptr) += (*sample_ptr++);
					(*prev_ptr++) /= 2.0;
				}

/* LINTED */
				errorStatus |= setMatrixValues(prevValues, toValueMatrix, toRow - 1, UNDEFINED_COUNT);

			}

			/*
			 * skip the next row
			 */

			toIndex_ptr[1] += 1.0;

			if ( (toIndex_ptr[1] + 1.0) == toAxisDimensions[1] )
			{
				/*
				 * set the last row in the plane equal
				 * to the second-to-the-last row
				 */

/* LINTED */
				errorStatus |= setMatrixValues(sampleValues, toValueMatrix, toRow + 1, UNDEFINED_COUNT);
			}

			prev_ptr = prevValues;
			prevValues = sampleValues;
			sampleValues = prev_ptr;

			toRow++;
		}

		/*
		 * increment the destination data array indices
		 * to move to the start of the next row
		 */

		for ( toAxis = 1; toAxis < toAxisCount; toAxis++ )
		{
			toIndex_ptr[toAxis] += 1.0;

			if ( toIndex_ptr[toAxis] < toAxisDimensions[toAxis] )
			{
				break;
			}

			/*
			 * reset this index to zero and
			 * increment the next index
			 */

			toIndex_ptr[toAxis] = 0.0;
		}

		toRow++;

	} /* end of while ( toRow < toRowCount ) */

  cleanupAndReturn:

	deallocateVector(&toIndexVector);
	deallocateMatrix(&indexTransform);
	deallocateVector(&fromIndexVector);
	(void)free(neighborOffsets);
	(void)free(neighborValues);
	(void)free(sampleValues);
	(void)free(prevValues);

	errorStatus |= setMatrixValuePtr((void *)NULL, toValueMatrix);
	deallocateMatrix(&toValueMatrix);

	if ( (errorStatus & cfiFailure) )
	{
		FAILURE_MSG("cannot map uniform quad-mesh values");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of mapUniQuadMeshValues() */



/************************************************************************
 *									*
 * getUniQuadMeshIndexTransform						*
 * ----------------------------						*
 *	Calculates the matrix that transforms floating-point		*
 *	indices	of one uniform quadrilateral mesh data array into	*
 *	floating-point indices of a second uniform quadrilateral	*
 *	mesh data array, such that the latter set of array indices	*
 *	corresponds to the same coordinate location as the former.	*
 *									*
 * NOTE-The coordinate axes of the first mesh MUST lie in the		*
 *	space spanned by the coordinate axes of the second mesh.	*
 *	This condition is NOT checked.					*
 *									*
 *	Assumes that the meshes have already been checked to have	*
 *	appropriate coordinate dimension and numbers of axes.		*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
getUniQuadMeshIndexTransform
(
	const cfiMesh *const meshA,		/* the first mesh            */
	const cfiMesh *const meshB,		/* the second mesh           */
	cfiMatrix     *const indexTransform	/* the transform of interest */
)

{
	cfiErrorStatus errorStatus;	/* the return value                  */

	cfiCounter n;			/* coordinate dimension of meshes    */
	cfiCounter i;			/* number of axes spanning meshA     */
	cfiCounter j;			/* number of axes spanning meshB     */

	cfiMatrix *Aprime;		/* see the extended comments below   */
	cfiMatrix *BprimeTranspose;
	cfiMatrix *L;
	cfiMatrix *Linverse;

	cfiVector *vector1;		/* coordinate vector storage         */
	cfiScalar *coord_ptr1;
	cfiVector *vector2;
	cfiScalar *coord_ptr2;
	cfiVector *vector3;
	cfiScalar *coord_ptr3;
	cfiVector *orthoVector;
	cfiScalar *ortho_ptr;

	cfiScalar epsilon;		/* a smallish number                 */

	cfiScalar scalar1;		/* intermediate value storage        */
	cfiScalar scalar2;

	cfiCounter axis;		/* loop control                      */
	cfiCounter orthoAxis;
	cfiCounter row;
	cfiCounter col;

	/*
	 * In the following,
	 *
	 *	x is an i-dimensional row vector of indices
	 *	  of the data array for meshA,
	 *
	 *	A is an ixn matrix containing the i n-dimensional
	 *	  linearly-independent coordinate axes for meshA,
	 *
	 *	a is an n-dimensional row vector containing
	 *	  the coordinate origin for meshA's vertices,
	 *
	 *	y is a j-dimensional row vector of indices
	 *	  of the data array for meshB,
	 *
	 *	B is a jxn matrix containing the j n-dimensional
	 *	  linearly-independent coordinate axes for meshB, and
	 *
	 *	b is an n-dimensional row vector containing
	 *	  the coordinate origin for meshB's vertices.
	 *
	 * If the coordinate axes of meshA lie in the space spanned by the
	 * coordinate axes of meshB, one can equate coordinate locations
	 * to obtain
	 *
	 *	 x A + a = y B + b.
	 *
	 * Rearranging and using homogeneous array indices for meshA,
	 * one has
	 *
	 *	x' A' = y B, where
	 *
	 * where
	 *
	 * 	x' = [ x | 1 ]
	 *
	 * and (page break follows)
	 *
	 *	     +-     -+
	 *	     |   A   |
	 *	A' = | ----- |
	 *	     |  a-b  |
	 *	     +-     -+.
	 *
	 * One can represent the matrix B as the product of a lower-
	 * triangular jxj matrix L and a jxn matrix B' containing
	 * j orthonormal coordinate axes spanning meshB, so that
	 *
	 *	x' A' = y L B'.
	 *
	 * Post-multiplying each side by the transpose of B' yields
	 *
	 *	        t
	 *	x' A' B'  = y L.
	 *
	 * The lower-triangular matrix L is inverted easily to yield
	 *
	 *	        t  -1
	 *	x' A' B'  L   = y.
	 *
	 * Thus, the (i+1)xj index transformation matrix of interest is
	 *
	 *	         t  -1
	 *	M = A' B'  L  .
	 */

	DEBUG_MSG(cfiTraceDebug, "getting uniform quad-mesh array index transform");

	n = getMeshCoordDimension(meshA);
	i = getMeshAxisCount(meshA);
	j = getMeshAxisCount(meshB);

	/*
	 * allocate workspace
	 */

	errorStatus = cfiSuccess;

	Aprime = (cfiMatrix *)NULL;
	BprimeTranspose = (cfiMatrix *)NULL;
	L = (cfiMatrix *)NULL;
	Linverse = (cfiMatrix *)NULL;
	vector1 = (cfiVector *)NULL;
	vector2 = (cfiVector *)NULL;
	vector3 = (cfiVector *)NULL;
	orthoVector = (cfiVector *)NULL;

	if ( ( (Aprime = allocateMatrix((cfiDataType)scalarType, i + 1, n))
		== (cfiMatrix *)NULL )
	||   ( (BprimeTranspose = allocateMatrix((cfiDataType)scalarType, j, n))
		== (cfiMatrix *)NULL )
	||   ( (L = allocateMatrix((cfiDataType)scalarType, j, j))
		== (cfiMatrix *)NULL )
	||   ( (Linverse = allocateMatrix((cfiDataType)scalarType, j, j))
		== (cfiMatrix *)NULL )
	||   ( (vector1 = allocateVector((cfiDataType)scalarType, n))
		== (cfiVector *)NULL )
	||   ( (coord_ptr1 = (cfiScalar *)getVectorValuePtr(vector1, 0))
		== (cfiScalar *)NULL )
	||   ( (vector2 = allocateVector((cfiDataType)scalarType, n))
		== (cfiVector *)NULL )
	||   ( (coord_ptr2 = (cfiScalar *)getVectorValuePtr(vector2, 0))
		== (cfiScalar *)NULL )
	||   ( (vector3 = allocateVector((cfiDataType)scalarType, n))
		== (cfiVector *)NULL )
	||   ( (coord_ptr3 = (cfiScalar *)getVectorValuePtr(vector3, 0))
		== (cfiScalar *)NULL )
	||   ( (orthoVector = allocateVector((cfiDataType)scalarType, n))
		== (cfiVector *)NULL )
	||   ( (ortho_ptr = (cfiScalar *)getVectorValuePtr(orthoVector, 0))
		== (cfiScalar *)NULL ) )
	{
		errorStatus |= cfiFailure;
		goto cleanupAndReturn;
	}

	/*
	 * calculate Aprime
	 */

	for ( axis = 0; axis < i; axis++ )
	{
		/*
		 * copy meshA axis to row in Aprime
		 */

		errorStatus |= getUniQuadMeshAxis(meshA, (cfiDataAxis)axis, coord_ptr1);
/* LINTED */
		errorStatus |= setMatrixValues(coord_ptr1, Aprime, axis, UNDEFINED_COUNT);
	}

	/*
	 * store translation from meshB vertex origin
	 * to meshA vertex origin in last row of Aprime
	 */

	errorStatus |= getUniQuadMeshOrigin(meshA, coord_ptr1);
	if ( getMeshType(meshA) == uniformQuadCenterType )
	{
		/*
		 * calculate position of [0][0]...[0]-th voxel center
		 */

		for ( axis = 0; axis < i; axis++ )
		{
			errorStatus |= getUniQuadMeshAxis(meshA, (cfiDataAxis)axis, coord_ptr3);
			errorStatus |= scaleVector(0.5, vector3);
			errorStatus |= vectorOpVector(vector1, (cfiOp)addOp, vector3, vector1);
		}
	}

	errorStatus |= getUniQuadMeshOrigin(meshB, coord_ptr2);
	if ( getMeshType(meshB) == uniformQuadCenterType )
	{
		/*
		 * calculate position of [0][0]...[0]-th voxel center
		 */

		for ( axis = 0; axis < j; axis++ )
		{
			errorStatus |= getUniQuadMeshAxis(meshB, (cfiDataAxis)axis, coord_ptr3);
			errorStatus |= scaleVector(0.5, vector3);
			errorStatus |= vectorOpVector(vector2, (cfiOp)addOp, vector3, vector2);
		}
	}

	errorStatus |= vectorOpVector(vector1, (cfiOp)subtractOp, vector2, vector1);
/* LINTED */
	errorStatus |= setMatrixValues(coord_ptr1, Aprime, i, UNDEFINED_COUNT);

	if ( (errorStatus & cfiFailure) )
	{
		goto cleanupAndReturn;
	}

	/*
	 * calculate L and BprimeTranspose
	 */

	epsilon = 0.000001;

	for ( axis = 0; axis < j; axis++ )
	{
		errorStatus |= getUniQuadMeshAxis(meshB, (cfiDataAxis)axis, coord_ptr1);
		errorStatus |= copyVector(vector1, vector2);

		for ( orthoAxis = 0; orthoAxis < axis; orthoAxis++ )
		{
			/*
			 * project meshB axis onto axis already orthonormalized
			 * and store projected component in L
			 */

/* LINTED */
			errorStatus |= getMatrixValues(BprimeTranspose, UNDEFINED_COUNT, orthoAxis, ortho_ptr);
			errorStatus |= vectorDotVector(vector1, orthoVector, &scalar1);
			errorStatus |= setMatrixValues(&scalar1, L, axis, orthoAxis);

			/*
			 * subtract projection from meshB axis
			 */

			errorStatus |= scaleVector(scalar1, orthoVector);
			errorStatus |= vectorOpVector(vector2, (cfiOp)subtractOp, orthoVector, vector2);

		} /* end of for ( orthoAxis = 0; ... ) loop */

		/*
		 * store orthonormal component in L and store
		 * orthonormalized axis in column of BprimeTranspose
		 */

		errorStatus |= vectorDotVector(vector2, vector2, &scalar1);
		if ( scalar1 < epsilon )
		{
			errorStatus |= cfiFailure;
			goto cleanupAndReturn;
		}
		scalar1 = sqrt(scalar1);
		errorStatus |= setMatrixValues(&scalar1, L, axis, axis);

		errorStatus |= scaleVector(1.0 / scalar1, vector2);
/* LINTED */
		errorStatus |= setMatrixValues(coord_ptr2, BprimeTranspose, UNDEFINED_COUNT, axis);

	} /* end of for ( axis = 0; ... ) loop */

	if ( (errorStatus & cfiFailure) )
	{
		goto cleanupAndReturn;
	}

	/*
	 * calculate Linverse
	 */

	for ( row = 0; row < j; row++ )
	{
/* LINTED */
		errorStatus |= getMatrixValues(L, row, UNDEFINED_COUNT, coord_ptr1);
		errorStatus |= getVectorValue(vector1, row, &scalar1);
		if ( fabs(scalar1) < epsilon )
		{
			errorStatus |= cfiFailure;
			goto cleanupAndReturn;
		}

		/*
		 * calculate lower diagonal elements in row
		 */

		for ( col = 0; col < row; col++ )
		{
/* LINTED */
			errorStatus |= getMatrixValues(Linverse, UNDEFINED_COUNT, col, coord_ptr2);
			errorStatus |= vectorDotVector(vector1, vector2, &scalar2);
			scalar2 /= (-scalar1);
			errorStatus |= setMatrixValues(&scalar2, Linverse, row, col);
		}

		/*
		 * calculate diagonal element in row
		 */

		scalar1 = 1.0 / scalar1;
		errorStatus |= setMatrixValues(&scalar1, Linverse, row, row);

	}

	if ( (errorStatus & cfiFailure) )
	{
		goto cleanupAndReturn;
	}

	/*
	 * calculate the index transform
	 */

	errorStatus |= matrixTimesMatrix(Aprime, BprimeTranspose, indexTransform);
	errorStatus |= matrixTimesMatrix(indexTransform, Linverse, indexTransform);

  cleanupAndReturn:

	deallocateMatrix(&Aprime);
	deallocateMatrix(&BprimeTranspose);
	deallocateMatrix(&L);
	deallocateMatrix(&Linverse);
	deallocateVector(&vector1);
	deallocateVector(&vector2);
	deallocateVector(&vector3);
	deallocateVector(&orthoVector);

	if ( (errorStatus & cfiFailure) )
	{
		FAILURE_MSG("cannot get mesh data array index transform");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of getUniQuadMeshIndexTransform() */



/************************************************************************
 *									*
 * sampleMeshValue							*
 * ---------------							*
 *	Samples a mesh data array at an arbitrary index location.	*
 *	Linearly interpolates between the values at the mesh vertices	*
 *	adjacent to the index location.					*
 *									*
 * NOTE-Does no error checking.						*
 *									*
 *	Assumes that indices have already been incremented by 0.5	*
 *	if (sub)nearestNeighbor sampling is to be done.			*
 *									*
 ********************************************************************bwr*/

void
sampleMeshValue
(
	const cfiMeshSamplingMethod samplingMethod,	/* sampling method    */

	const cfiCounter   axisCount,		/* # of axes spanning mesh    */
	const cfiCounter  *const axisDimensions,/* dimensions of axes         */
	const cfiScalar   *const scalarIndices,	/* indices of sample location */
	const cfiDataType  valueType,		/* type of mesh values        */
	void              *value_ptr,		/* pointer to mesh values     */
	cfiCounter         neighborCount,	/* # of corners in mesh voxel */
	const cfiCounter  *neighborOffsets,	/* pointer offsets to corners */
	cfiScalar         *const neighborValues,/* voxel corner value storage */
	cfiScalar         *const sampleValue	/* sampled value of interest  */
)

{
        cfiScalar scalarIndex;		/* floating-point index for location  */
	cfiScalar u;			/* linear interpolation weight factor */

	cfiCounter neighborZeroOffset;	/* pointer offset to zeroth neighbor  */

	cfiScalar n00, n01, n02, n03;
	cfiScalar n04, n05, n06, n07;
	cfiScalar n08, n09, n10, n11;
	cfiScalar n12, n13, n14, n15;

	cfiScalar  *interp_ptr;		/* loop control                       */
	cfiScalar  *neighbor_ptr;
	cfiScalar  *end_neighbor;
	cfiCounter  axis;

	/*
	 * check the floating point indices and calculate
	 * the pointer offset to the zero-th neighbor
	 */

	neighborZeroOffset = 0;
	for ( axis = 0; axis < axisCount; axis++ )
	{
		scalarIndex = scalarIndices[axis];

		if ( ( scalarIndex <= 0.0 )
		||   ( scalarIndex >= (axisDimensions[axis] - 1) ) )
		{
			/*
			 * the floating-point index along this axis
			 * takes us outside of the coordinate domain
			 * of the mesh whose data array is to be sampled
			 */

			(*sampleValue) = 0.0;
			return;
		}

		/*
		 * truncate the index and update the offset
		 */

		neighborZeroOffset += ((cfiCounter)scalarIndex) * neighborOffsets[1L<<axis];
	}

	if ( samplingMethod != (cfiMeshSamplingMethod)linearInterp )
	{
		switch (valueType)
		{
		  case byteType:

			(*sampleValue) = (cfiScalar)(*(((unsigned char *)value_ptr) + neighborZeroOffset));
			return;

		  case shortType:

			(*sampleValue) = (cfiScalar)(*(((int *)value_ptr) + neighborZeroOffset));
			return;

		  case longType:

			(*sampleValue) = (cfiScalar)(*(((long *)value_ptr) + neighborZeroOffset));
			return;

		  case floatType:

			(*sampleValue) = (cfiScalar)(*(((float *)value_ptr) + neighborZeroOffset));
			return;

		  case doubleType:

			(*sampleValue) = (cfiScalar)(*(((double *)value_ptr) + neighborZeroOffset));
			return;

		  default:

			(*sampleValue) = 0.0;
			return;

		} /* end of switch (valueType) */

	} /* end of if ( samplingMethod != linearInterp ) statement */

	/*
	 * copy the neighboring values to the voxel corner value storage
	 */

	neighbor_ptr = neighborValues;
	end_neighbor = neighborValues + neighborCount;

	switch (valueType)
	{
	  case byteType:

		value_ptr = (void *)(((unsigned char *)value_ptr) + neighborZeroOffset);

		switch (axisCount)
		{
		  case 4:

			n15 = (cfiScalar)(*(((unsigned char *)value_ptr) + neighborOffsets[15]));
			n14 = (cfiScalar)(*(((unsigned char *)value_ptr) + neighborOffsets[14]));
			n13 = (cfiScalar)(*(((unsigned char *)value_ptr) + neighborOffsets[13]));
			n12 = (cfiScalar)(*(((unsigned char *)value_ptr) + neighborOffsets[12]));
			n11 = (cfiScalar)(*(((unsigned char *)value_ptr) + neighborOffsets[11]));
			n10 = (cfiScalar)(*(((unsigned char *)value_ptr) + neighborOffsets[10]));
			n09 = (cfiScalar)(*(((unsigned char *)value_ptr) + neighborOffsets[9]));
			n08 = (cfiScalar)(*(((unsigned char *)value_ptr) + neighborOffsets[8]));

/* FALLTHROUGH */
		  case 3:

			n07 = (cfiScalar)(*(((unsigned char *)value_ptr) + neighborOffsets[7]));
			n06 = (cfiScalar)(*(((unsigned char *)value_ptr) + neighborOffsets[6]));
			n05 = (cfiScalar)(*(((unsigned char *)value_ptr) + neighborOffsets[5]));
			n04 = (cfiScalar)(*(((unsigned char *)value_ptr) + neighborOffsets[4]));

/* FALLTHROUGH */
		  case 2:

			n03 = (cfiScalar)(*(((unsigned char *)value_ptr) + neighborOffsets[3]));
			n02 = (cfiScalar)(*(((unsigned char *)value_ptr) + neighborOffsets[2]));

/* FALLTHROUGH */
		  case 1:

			n01 = (cfiScalar)(*(((unsigned char *)value_ptr) + neighborOffsets[1]));
			n00 = (cfiScalar)(*(((unsigned char *)value_ptr) + neighborOffsets[0]));

			break;

		  default:

			while ( neighbor_ptr < end_neighbor )
			{
				(*neighbor_ptr++) = (cfiScalar)(*(((unsigned char *)value_ptr) + (*neighborOffsets++)));
			}

			break;

		} /* end of switch (axisCount) */

		break;

	  case shortType:

		value_ptr = (void *)(((int *)value_ptr) + neighborZeroOffset);

		switch (axisCount)
		{
		  case 4:

			n15 = (cfiScalar)(*(((int *)value_ptr) + neighborOffsets[15]));
			n14 = (cfiScalar)(*(((int *)value_ptr) + neighborOffsets[14]));
			n13 = (cfiScalar)(*(((int *)value_ptr) + neighborOffsets[13]));
			n12 = (cfiScalar)(*(((int *)value_ptr) + neighborOffsets[12]));
			n11 = (cfiScalar)(*(((int *)value_ptr) + neighborOffsets[11]));
			n10 = (cfiScalar)(*(((int *)value_ptr) + neighborOffsets[10]));
			n09 = (cfiScalar)(*(((int *)value_ptr) + neighborOffsets[9]));
			n08 = (cfiScalar)(*(((int *)value_ptr) + neighborOffsets[8]));

/* FALLTHROUGH */
		  case 3:

			n07 = (cfiScalar)(*(((int *)value_ptr) + neighborOffsets[7]));
			n06 = (cfiScalar)(*(((int *)value_ptr) + neighborOffsets[6]));
			n05 = (cfiScalar)(*(((int *)value_ptr) + neighborOffsets[5]));
			n04 = (cfiScalar)(*(((int *)value_ptr) + neighborOffsets[4]));

/* FALLTHROUGH */
		  case 2:

			n03 = (cfiScalar)(*(((int *)value_ptr) + neighborOffsets[3]));
			n02 = (cfiScalar)(*(((int *)value_ptr) + neighborOffsets[2]));

/* FALLTHROUGH */
		  case 1:

			n01 = (cfiScalar)(*(((int *)value_ptr) + neighborOffsets[1]));
			n00 = (cfiScalar)(*(((int *)value_ptr) + neighborOffsets[0]));

			break;

		  default:

			while ( neighbor_ptr < end_neighbor )
			{
				(*neighbor_ptr++) = (cfiScalar)(*(((int *)value_ptr) + (*neighborOffsets++)));
			}

			break;

		} /* end of switch (axisCount) */

		break;

	  case longType:

		value_ptr = (void *)(((long *)value_ptr) + neighborZeroOffset);

		switch (axisCount)
		{
		  case 4:

			n15 = (cfiScalar)(*(((long *)value_ptr) + neighborOffsets[15]));
			n14 = (cfiScalar)(*(((long *)value_ptr) + neighborOffsets[14]));
			n13 = (cfiScalar)(*(((long *)value_ptr) + neighborOffsets[13]));
			n12 = (cfiScalar)(*(((long *)value_ptr) + neighborOffsets[12]));
			n11 = (cfiScalar)(*(((long *)value_ptr) + neighborOffsets[11]));
			n10 = (cfiScalar)(*(((long *)value_ptr) + neighborOffsets[10]));
			n09 = (cfiScalar)(*(((long *)value_ptr) + neighborOffsets[9]));
			n08 = (cfiScalar)(*(((long *)value_ptr) + neighborOffsets[8]));

/* FALLTHROUGH */
		  case 3:

			n07 = (cfiScalar)(*(((long *)value_ptr) + neighborOffsets[7]));
			n06 = (cfiScalar)(*(((long *)value_ptr) + neighborOffsets[6]));
			n05 = (cfiScalar)(*(((long *)value_ptr) + neighborOffsets[5]));
			n04 = (cfiScalar)(*(((long *)value_ptr) + neighborOffsets[4]));

/* FALLTHROUGH */
		  case 2:

			n03 = (cfiScalar)(*(((long *)value_ptr) + neighborOffsets[3]));
			n02 = (cfiScalar)(*(((long *)value_ptr) + neighborOffsets[2]));

/* FALLTHROUGH */
		  case 1:

			n01 = (cfiScalar)(*(((long *)value_ptr) + neighborOffsets[1]));
			n00 = (cfiScalar)(*(((long *)value_ptr) + neighborOffsets[0]));

			break;

		  default:

			while ( neighbor_ptr < end_neighbor )
			{
				(*neighbor_ptr++) = (cfiScalar)(*(((long *)value_ptr) + (*neighborOffsets++)));
			}

			break;

		} /* end of switch (axisCount) */

		break;

	  case floatType:

		value_ptr = (void *)(((float *)value_ptr) + neighborZeroOffset);

		switch (axisCount)
		{
		  case 4:

			n15 = (cfiScalar)(*(((float *)value_ptr) + neighborOffsets[15]));
			n14 = (cfiScalar)(*(((float *)value_ptr) + neighborOffsets[14]));
			n13 = (cfiScalar)(*(((float *)value_ptr) + neighborOffsets[13]));
			n12 = (cfiScalar)(*(((float *)value_ptr) + neighborOffsets[12]));
			n11 = (cfiScalar)(*(((float *)value_ptr) + neighborOffsets[11]));
			n10 = (cfiScalar)(*(((float *)value_ptr) + neighborOffsets[10]));
			n09 = (cfiScalar)(*(((float *)value_ptr) + neighborOffsets[9]));
			n08 = (cfiScalar)(*(((float *)value_ptr) + neighborOffsets[8]));

/* FALLTHROUGH */
		  case 3:

			n07 = (cfiScalar)(*(((float *)value_ptr) + neighborOffsets[7]));
			n06 = (cfiScalar)(*(((float *)value_ptr) + neighborOffsets[6]));
			n05 = (cfiScalar)(*(((float *)value_ptr) + neighborOffsets[5]));
			n04 = (cfiScalar)(*(((float *)value_ptr) + neighborOffsets[4]));

/* FALLTHROUGH */
		  case 2:

			n03 = (cfiScalar)(*(((float *)value_ptr) + neighborOffsets[3]));
			n02 = (cfiScalar)(*(((float *)value_ptr) + neighborOffsets[2]));

/* FALLTHROUGH */
		  case 1:

			n01 = (cfiScalar)(*(((float *)value_ptr) + neighborOffsets[1]));
			n00 = (cfiScalar)(*(((float *)value_ptr) + neighborOffsets[0]));

			break;

		  default:

			while ( neighbor_ptr < end_neighbor )
			{
				(*neighbor_ptr++) = (cfiScalar)(*(((float *)value_ptr) + (*neighborOffsets++)));
			}

			break;

		} /* end of switch (axisCount) */

		break;

	  case doubleType:

		value_ptr = (void *)(((double *)value_ptr) + neighborZeroOffset);

		switch (axisCount)
		{
		  case 4:

			n15 = (cfiScalar)(*(((double *)value_ptr) + neighborOffsets[15]));
			n14 = (cfiScalar)(*(((double *)value_ptr) + neighborOffsets[14]));
			n13 = (cfiScalar)(*(((double *)value_ptr) + neighborOffsets[13]));
			n12 = (cfiScalar)(*(((double *)value_ptr) + neighborOffsets[12]));
			n11 = (cfiScalar)(*(((double *)value_ptr) + neighborOffsets[11]));
			n10 = (cfiScalar)(*(((double *)value_ptr) + neighborOffsets[10]));
			n09 = (cfiScalar)(*(((double *)value_ptr) + neighborOffsets[9]));
			n08 = (cfiScalar)(*(((double *)value_ptr) + neighborOffsets[8]));

/* FALLTHROUGH */
		  case 3:

			n07 = (cfiScalar)(*(((double *)value_ptr) + neighborOffsets[7]));
			n06 = (cfiScalar)(*(((double *)value_ptr) + neighborOffsets[6]));
			n05 = (cfiScalar)(*(((double *)value_ptr) + neighborOffsets[5]));
			n04 = (cfiScalar)(*(((double *)value_ptr) + neighborOffsets[4]));

/* FALLTHROUGH */
		  case 2:

			n03 = (cfiScalar)(*(((double *)value_ptr) + neighborOffsets[3]));
			n02 = (cfiScalar)(*(((double *)value_ptr) + neighborOffsets[2]));

/* FALLTHROUGH */
		  case 1:

			n01 = (cfiScalar)(*(((double *)value_ptr) + neighborOffsets[1]));
			n00 = (cfiScalar)(*(((double *)value_ptr) + neighborOffsets[0]));

			break;

		  default:

			while ( neighbor_ptr < end_neighbor )
			{
				(*neighbor_ptr++) = (cfiScalar)(*(((double *)value_ptr) + (*neighborOffsets++)));
			}

			break;

		} /* end of switch (axisCount) */

		break;

	  default:

		(*sampleValue) = 0.0;
		return;

	} /* end of switch (valueType) */

	/*
	 * (page break follows)
	 * linearly interpolate between neighboring values along each axis
	 * (e.g. in 3-D one has
	 *
	 *		      n6          n3'         n7
	 *			+-----------*-----------+
	 *		                   /
	 *		      /           /           /
	 *		        |        /              |
	 *		    /           /           /
	 *		        |      * n1''           |
	 *		  /           /|          /
	 *	      n4        | n2'/ |      n5        |
	 *		+-----------*-----------+
	 *		        |      |                |
	 *		|              |        |
	 *		        +------|----*-----------+
	 *		|        n2    |   / n1'|        n3
	 *		      /      S *  /           /
	 *		|              | /      |
	 *		    /          |/           /
	 *		|              * n0''   |
	 *		  /           /           /
	 *                           /
	 *		+-----------*-----------+
	 *               n0          n0'         n1
	 *
	 *	n0' = n0 + (u * (n1 - n0));
	 *	n1' = n2 + (u * (n3 - n2));
	 *	n2' = n4 + (u * (n5 - n4));
	 *	n3' = n6 + (u * (n7 - n6));
	 *
	 *	n0'' = n0' + (v * (n1' - n0'));
	 *	n1'' = n2' + (v * (n3' - n2'));
	 *
	 *	S = n0'' + (w * (n1'' - n0''));
	 *
	 * where [u, v, w] are the fractional parts of the three
	 * floating-point indices for the coordinate location to be
	 * sampled; n0, n1, n2, n3, n4, n5, n6, and n7 are the mesh
	 * data array values at the eight neighboring vertices; and
	 * S is the linearly-interpolated sample value)
	 */

	switch (axisCount)
	{

	  /*
	   * NOTE--The default code below processes indices in the
	   *       order u,v,w,... as described above.  The optimized
	   *       special cases process indices "backwards" in the
	   *       order ...,w,v,u.
	   */

	  case 4:

		u = scalarIndices[3];
		u -= (cfiCounter)u;

		n00 += u * (n08 - n00);
		n01 += u * (n09 - n01);
		n02 += u * (n10 - n02);
		n03 += u * (n11 - n03);
		n04 += u * (n12 - n04);
		n05 += u * (n13 - n05);
		n06 += u * (n14 - n06);
		n07 += u * (n15 - n07);

/* FALLTHROUGH */
	  case 3:

		u = scalarIndices[2];
		u -= (cfiCounter)u;

		n00 += u * (n04 - n00);
		n01 += u * (n05 - n01);
		n02 += u * (n06 - n02);
		n03 += u * (n07 - n03);

/* FALLTHROUGH */
	  case 2:

		u = scalarIndices[1];
		u -= (cfiCounter)u;

		n00 += u * (n02 - n00);
		n01 += u * (n03 - n01);

/* FALLTHROUGH */
	  case 1:

		u = scalarIndices[0];
		u -= (cfiCounter)u;

		(*sampleValue) = n00 + (u * (n01 - n00));
		return;

	  default:

		for ( axis = 0; axis < axisCount; axis++ )
		{
			u = scalarIndices[axis];
			u -= (cfiCounter)u;

			interp_ptr = neighborValues;
			neighbor_ptr = neighborValues;
			end_neighbor = neighborValues + neighborCount;

			while ( neighbor_ptr < end_neighbor )
			{
				(*interp_ptr++) = (*neighbor_ptr) + (u * ((*(neighbor_ptr + 1)) - (*neighbor_ptr)));
				neighbor_ptr += 2;
			}

			neighborCount /= 2;
		}

		(*sampleValue) = (*neighborValues);
		return;

	} /* end of switch (axisCount) */

} /* end of sampleMeshValue() */



/************************************************************************
 *									*
 * allocateMesh								*
 * ------------								*
 *	Returns a pointer to a mesh of the specified geometry		*
 *	type, vertex coordinate dimension, and data value type		*
 *	and dimensions.  The vertex coordinates and data values		*
 *	are set to zero.						*
 *									*
 * NOTE-Passing a NULL "dimensions" pointer is equivalent		*
 *	to specifying that there are zero vertices along		*
 *	each mesh axis.							*
 *									*
 ********************************************************************bwr*/

cfiMesh *
allocateMesh
(
	const cfiMeshType  meshType,		/* type of mesh geometry      */
	const cfiCounter   coordCount,		/* coordinate dimension       */
	const cfiDataType  dataType,		/* type of data values        */
	const cfiCounter   axisCount,		/* # of axes spanning mesh    */
	const cfiCounter  *const dimensions	/* #'s of vertices along axes */
)

{
	cfiMesh *mesh;				/* the return value           */

	DEBUG_MSG(cfiMemoryDebug, "allocating a mesh");

	if ( (    ( meshType != (cfiMeshType)uniformQuadCornerType )
	       && ( meshType != (cfiMeshType)uniformQuadCenterType ) )
	||   ( axisCount <= 0 )
	||   ( axisCount > coordCount )
/* LINTED */
	||   ( (mesh = (cfiMesh *)malloc(sizeof(cfiMesh)))
		== (cfiMesh *)NULL ) )
	{
		FAILURE_MSG("cannot allocate mesh");
		return ((cfiMesh *)NULL);
	}

	mesh->meshType = meshType;
	mesh->meshCoords = (cfiMatrix *)NULL;
	mesh->meshValues = (cfiArray *)NULL;

	if ( ( (mesh->meshCoords = allocateMatrix((cfiDataType)scalarType, axisCount + 1, coordCount))
		== (cfiMatrix *)NULL )
	||   ( (mesh->meshValues = allocateArray(dataType, axisCount, dimensions))
		== (cfiArray *)NULL ) )
	{
		FAILURE_MSG("cannot allocate mesh");
		deallocateMesh(&mesh);
		return ((cfiMesh *)NULL);
	}

	return (mesh);

} /* end of allocateMesh() */



/************************************************************************
 *									*
 * deallocateMesh							*
 * --------------							*
 *	Deallocates memory for a mesh.					*
 *									*
 ********************************************************************bwr*/

void
deallocateMesh
(
	cfiMesh **const mesh_ptr	/* the mesh of interest */
)

{
	DEBUG_MSG(cfiMemoryDebug, "deallocating mesh");

	if ( (*mesh_ptr) == (cfiMesh *)NULL )
	{
		return;
	}

	deallocateMatrix(&((*mesh_ptr)->meshCoords));
	deallocateArray(&((*mesh_ptr)->meshValues));
	(void)free(*mesh_ptr);
	(*mesh_ptr) = (cfiMesh *)NULL;

} /* end of deallocateMesh() */



/************************************************************************
 *									*
 * writeMesh								*
 * ---------								*
 *	Writes a mesh to a stream.					*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
writeMesh
(
	const cfiMesh *const mesh,	/* the mesh of interest */
	FILE          *const stream	/* the output stream    */
)

{
	cfiErrorStatus errorStatus;	/* the return value     */

	DEBUG_MSG(cfiInOutDebug, "writing mesh");

	if ( mesh == (cfiMesh *)NULL )
	{
		if ( fprintf(stream, "null mesh\n") < 0 )
		{
			FAILURE_MSG("cannot write mesh");
			return (cfiFailure);
		}

		return (cfiSuccess);
	}

	switch (getMeshType(mesh))
	{
	  case uniformQuadCornerType:

		(void)fprintf(stream, "uniformQuadCornerType mesh\n");
		break;

	  case uniformQuadCenterType:

		(void)fprintf(stream, "uniformQuadCenterType mesh\n");
		break;

	  default:

		(void)fprintf(stream, "unknownMeshType\n");
		break;
	}

	errorStatus = cfiSuccess;

	(void)fprintf(stream, "coordinate ");
	errorStatus |= writeMatrix(mesh->meshCoords, stream);

	(void)fprintf(stream, "data value ");
	errorStatus |= writeArray(mesh->meshValues, stream);

	if ( (errorStatus & cfiFailure) )
	{
		FAILURE_MSG("cannot write mesh");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of writeMesh() */



/************************************************************************
 *									*
 * checkMesh								*
 * ---------								*
 *	Checks if a mesh has data values at				*
 *	and coordinates for its vertices.				*
 *									*
 ********************************************************************bwr*/

cfiMeshStatus
checkMesh
(
	const cfiMesh *const mesh	/* the mesh of interest             */
)

{
	cfiMeshStatus meshStatus;	/* the return value                 */

	cfiMeshType     meshType;	/* the type of the mesh             */
	cfiArrayStatus  valueStatus;	/* status of the data values        */
	cfiMatrixStatus coordStatus;	/* status of the coordinates        */

	long axisCount;			/* number of axes spanning the mesh */

	DEBUG_MSG(cfiTraceDebug, "checking mesh");

	meshStatus = emptyMesh;

	if ( (meshType = getMeshType(mesh))
		== unknownMeshType )
	{
		return (meshStatus);
	}

	/*
	 * check data values
	 */

	valueStatus = checkArray(mesh->meshValues);

	if ( (valueStatus & arrayHasValues) )
	{
		meshStatus |= meshHasValues;
	}

	if ( (valueStatus & arrayHasDimensions) )
	{
		meshStatus |= meshHasDimensions;
	}
	else
	{
		return (meshStatus);
	}

	/*
	 * check coordinates
	 */

	coordStatus = checkMatrix(mesh->meshCoords);

	if ( ( ! (coordStatus & matrixHasValues) )
	||   ( ! (coordStatus & matrixHasDimensions) )
	||   ( getMatrixValueType(mesh->meshCoords) != scalarType )
	||   ( (axisCount = getMeshAxisCount(mesh))
		<= 0 )
	||   ( axisCount > getMeshCoordDimension(mesh) ) )
	{
		return (meshStatus);
	}

	if ( (    ( meshType == uniformQuadCornerType )
	       || ( meshType == uniformQuadCenterType ) )
	&&   ( getMatrixDimension(mesh->meshCoords, (cfiDataAxis)rowAxis) == (axisCount + 1) ) )
	{
		meshStatus |= meshHasCoords;
	}

	return (meshStatus);

} /* end of checkMesh() */



/************************************************************************
 *									*
 * getUniQuadMeshOrigin							*
 * --------------------							*
 *	Returns the position vector for the origin of a uniform		*
 *	quadrilateral mesh, in a contiguous block of scalars.		*
 *									*
 * NOTE-For a uniformQuadCornerType mesh, the [0][0]...[0]-th		*
 *	mesh vertex is located at the origin.  Otherwise, for		*
 *	a uniformQuadCenterType mesh, the [0][0]...[0]-th mesh		*
 *	vertex is displaced from the origin in the positive		*
 *	direction along each axis by an amount equal to one-half	*
 *	of the vertex spacing for the axis.				*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
getUniQuadMeshOrigin
(
	const cfiMesh *const mesh,	/* the mesh of interest             */
	cfiScalar     *coord_ptr	/* destination of the vector coords */
)

{
	cfiCounter axisCount;		/* number of axes spanning the mesh */

	/*
	 * the mesh origin is stored in the last row
	 * of the coordinate matrix (the preceding rows
	 * store the mesh axis directions)
	 */

	if ( ( (axisCount = getMeshAxisCount(mesh))
		<= 0 )
	||   (    ( mesh->meshType != uniformQuadCornerType )
	       && ( mesh->meshType != uniformQuadCenterType ) )
/* LINTED */
	||   ( getMatrixValues(mesh->meshCoords, axisCount, UNDEFINED_COUNT, coord_ptr) != cfiSuccess ) )
	{
		FAILURE_MSG("cannot get uniform quad-mesh origin");
		return (cfiFailure);
	}

	return (cfiSuccess);

} /* end of getUniQuadMeshOrigin() */



/************************************************************************
 *									*
 * setUniQuadMeshOrigin							*
 * --------------------							*
 *	Sets the position vector for the origin of a uniform		*
 *	quadrilateral mesh, from a contiguous block of scalars.		*
 *									*
 * NOTE-For a uniformQuadCornerType mesh, the origin should be		*
 *	located at [0][0]...[0]-th mesh vertex.  Otherwise, for		*
 *	a uniformQuadCenterType mesh, the origin should be		*
 *	displaced from the [0][0]...[0]-th mesh vertex in the		*
 *	negative direction along each axis by an amount equal		*
 *	to one-half of the vertex spacing for that axis.		*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
setUniQuadMeshOrigin
(
	const cfiScalar *coord_ptr,	/* source of the vector coords      */
	cfiMesh         *const mesh	/* the mesh of interest             */
)

{
	cfiCounter axisCount;		/* number of axes spanning the mesh */

	/*
	 * the mesh origin is stored in the last row
	 * of the coordinate matrix (the preceding rows
	 * store the mesh axis directions)
	 */

	if ( ( (axisCount = getMeshAxisCount(mesh))
		<= 0 )
	||   (    ( mesh->meshType != uniformQuadCornerType )
	       && ( mesh->meshType != uniformQuadCenterType ) )
/* LINTED */
	||   ( setMatrixValues(coord_ptr, mesh->meshCoords, axisCount, UNDEFINED_COUNT) != cfiSuccess ) )
	{
		FAILURE_MSG("cannot set uniform quad-mesh origin");
		return (cfiFailure);
	}

	return (cfiSuccess);

} /* end of setUniQuadMeshOrigin() */



/************************************************************************
 *									*
 * getUniQuadMeshAxis							*
 * ------------------							*
 *	Returns the direction vector for an axis of a uniform		*
 *	quadrilateral mesh, in a contiguous block of scalars.		*
 *									*
 * NOTE-The magnitude of the direction vector is equal			*
 *	to the spacing between vertices along the axis.			*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
getUniQuadMeshAxis
(
	const cfiMesh     *const mesh,	/* the mesh of interest             */
	const cfiDataAxis  axis,	/* the axis of interest             */
	cfiScalar         *coord_ptr	/* destination of the vector coords */
)

{
	/*
	 * the mesh axis directions are stored in
	 * the first rows of the coordinate matrix
	 * (the last row stores the mesh origin)
	 */

	if ( ( axis < 0 )
	||   ( axis >= getMeshAxisCount(mesh) )
	||   (    ( mesh->meshType != uniformQuadCornerType )
	       && ( mesh->meshType != uniformQuadCenterType ) )
/* LINTED */
	||   ( getMatrixValues(mesh->meshCoords, axis, UNDEFINED_COUNT, coord_ptr) != cfiSuccess ) )
	{
		FAILURE_MSG("cannot get uniform quad-mesh axis direction");
		return (cfiFailure);
	}

	return (cfiSuccess);

} /* end of getUniQuadMeshAxis() */



/************************************************************************
 *									*
 * setUniQuadMeshAxis							*
 * ------------------							*
 *	Sets the direction vector for an axis of a uniform		*
 *	quadrilateral mesh, from a contiguous block of scalars.		*
 *									*
 * NOTE-The magnitude of the direction vector should be equal		*
 *	to the spacing between vertices along the axis.			*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
setUniQuadMeshAxis
(
	const cfiScalar   *coord_ptr,	/* source of the vector coords */
	cfiMesh           *const mesh,	/* the mesh of interest        */
	const cfiDataAxis  axis		/* the axis of interest        */
)

{
	/*
	 * the mesh axis directions are stored in
	 * the first rows of the coordinate matrix
	 * (the last row stores the mesh origin)
	 */

	if ( ( axis < 0 )
	||   ( axis >= getMeshAxisCount(mesh) )
	||   (    ( mesh->meshType != uniformQuadCornerType )
	       && ( mesh->meshType != uniformQuadCenterType ) )
/* LINTED */
	||   ( setMatrixValues(coord_ptr, mesh->meshCoords, axis, UNDEFINED_COUNT) != cfiSuccess ) )
	{
		FAILURE_MSG("cannot set uniform quad-mesh axis direction");
		return (cfiFailure);
	}

	return (cfiSuccess);

} /* end of setUniQuadMeshAxis() */



/************************************************************************
 *									*
 * getMeshType								*
 * -----------								*
 *	Returns the type of geometry for a mesh.			*
 *									*
 ********************************************************************bwr*/

cfiMeshType
getMeshType
(
	const cfiMesh *const mesh	/* the mesh of interest */
)

{
	cfiMeshType meshType;		/* the return value     */

	if ( mesh == (cfiMesh *)NULL )
	{
		DEBUG_MSG(cfiReturnDebug, "type of mesh is unknown");
		return (unknownMeshType);
	}

	meshType = mesh->meshType;

	if ( ( meshType == uniformQuadCornerType )
	||   ( meshType == uniformQuadCenterType ) )
	{
		return (meshType);
	}

	DEBUG_MSG(cfiReturnDebug, "type of mesh is unknown");
	return (unknownMeshType);

} /* end of getMeshType() */



/************************************************************************
 *									*
 * getMeshCoordDimension						*
 * ---------------------						*
 *	Returns the vertex coordinate dimension for a mesh.		*
 *									*
 ********************************************************************bwr*/

cfiCounter
getMeshCoordDimension
(
	const cfiMesh *const mesh	/* the mesh of interest */
)

{
	cfiCounter coordDimension;	/* the return value     */

	if ( ( mesh == (cfiMesh *)NULL )
	||   ( (coordDimension = getMatrixDimension(mesh->meshCoords, (cfiDataAxis)columnAxis))
		<= 0 ) )
	{
		DEBUG_MSG(cfiReturnDebug, "mesh coordinate dimension is ill-defined");
		return (0);
	}

	return (coordDimension);

} /* end of getMeshCoordDimension() */



/************************************************************************
 *									*
 * getMeshAxisCount							*
 * ----------------							*
 *	Returns the number of axes spanning a mesh.			*
 *									*
 ********************************************************************bwr*/

cfiCounter
getMeshAxisCount
(
	const cfiMesh *const mesh	/* the mesh of interest */
)

{
	cfiCounter axisCount;		/* the return value     */

	if ( ( mesh == (cfiMesh *)NULL )
	||   ( (axisCount = getArrayAxisCount(mesh->meshValues))
		<= 0 ) )
	{
		DEBUG_MSG(cfiReturnDebug, "number of mesh axes is ill-defined");
		return (0);
	}

	return (axisCount);

} /* end of getMeshAxisCount() */



/************************************************************************
 *									*
 * getMeshDataArrayPtr							*
 * -------------------							*
 *	Returns a pointer to the data array for a mesh.			*
 *									*
 ********************************************************************bwr*/

cfiArray *
getMeshDataArrayPtr
(
	const cfiMesh *const mesh	/* the mesh of interest */
)

{
	if ( ( mesh == (cfiMesh *)NULL )
	||   ( mesh->meshValues == (cfiArray *)NULL ) )
	{
		DEBUG_MSG(cfiReturnDebug, "cannot get pointer to mesh data array");
		return ((cfiArray *)NULL);
	}

	return (mesh->meshValues);

} /* end of getMeshDataArrayPtr() */
