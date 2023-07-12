/*
 * SCCS:  @(#)cfiArray.c  1.5  6/7/95  11:33:47
 */

/************************************************************************
 *									*
 * -----------------							*
 * FUNCTIONS HEREIN:							*
 * -----------------							*
 *									*
 * getArrayValues							*
 * --------------							*
 *	Retrieves values from an array, and stores the values		*
 *	in a contiguous block of scalars.				*
 *									*
 * setArrayValues							*
 * --------------							*
 *	Converts a contiguous block of scalars to the type of		*
 *	an array, and stores the resulting values in the array.		*
 *									*
 * allocateArray							*
 * -------------							*
 *	Returns a pointer to an array of the specified type		*
 *	and dimensions, with its values set to zero.			*
 *									*
 * deallocateArray							*
 * ---------------							*
 *	Deallocates memory for an array.				*
 *									*
 * arrayOpArray								*
 * ------------								*
 *	Performs pointwise arithmetic operations on two arrays.		*
 *									*
 * scaleArray								*
 * ----------								*
 *	Multiplies an array by a scalar.				*
 *									*
 * copyArray								*
 * ---------								*
 *	Copies values from an array of arbitrary type			*
 *	to another array of arbitrary type.				*
 *									*
 * writeArray								*
 * ----------								*
 *	Writes an array to a stream.					*
 *									*
 * checkArray								*
 * ----------								*
 *	Checks if an array has values of a known type			*
 *	and valid dimensions.						*
 *									*
 * zeroArray								*
 * ---------								*
 *	Sets the values in an array to zero.				*
 *									*
 * getArrayValueType							*
 * -----------------							*
 *	Returns the type of values in an array.				*
 *									*
 * setArrayValueType							*
 * -----------------							*
 *	Sets the type of values in an array.				*
 *									*
 * getArrayAxisCount							*
 * -----------------							*
 *	Returns the number of axes for an array.			*
 *									*
 * setArrayAxisCount							*
 * -----------------							*
 *	Sets the number of axes for an array.				*
 *									*
 * getArrayDimension							*
 * -----------------							*
 *	Returns the number of values along an array axis.		*
 *									*
 * setArrayDimension							*
 * -----------------							*
 *	Sets the number of values along an array axis.			*
 *									*
 * getArrayValuePtr							*
 * ----------------							*
 *	Returns a void pointer to a value in an array.			*
 *									*
 * setArrayValuePtr							*
 * ----------------							*
 *	Sets the pointer to the values in an array.			*
 *									*
 ************************************************************************/



#include <stdio.h>

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



static char l_msg[MSG_LENGTH];		/* for compound messages */



/************************************************************************
 *									*
 * getArrayValues							*
 * --------------							*
 *	Retrieves values from an array, and stores the values		*
 *	in a contiguous block of scalars.				*
 *									*
 * NOTE-Passing a contiguous set of UNDEFINED_COUNT values		*
 *	in "indices" gets all of the values from a subarray.		*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
getArrayValues
(
	const cfiArray   *const array,		/* the source array          */
	const cfiCounter *const indices,	/* array indices of interest */
	cfiScalar        *value_ptr		/* destination of the values */
)

{
	cfiErrorStatus errorStatus;	/* the return value                  */

	cfiArrayStatus arrayStatus;	/* status of the array               */

	cfiCounter axisCount;		/* number of array axes              */
	cfiCounter subStart[MAX_AXIS];	/* indices for start of subarray     */

	cfiMatrix  *subMatrix;		/* simplifies access to subarray     */
	cfiCounter  rowCount;		/* number of rows in submatrix       */
	cfiCounter  colCount;		/* number of columns in submatrix    */
	cfiCounter  i;			/* loop control                      */

	if ( (g_debugStatus & cfiArgDebug) )
	{
		arrayStatus = checkArray(array);

		if ( ( ! (arrayStatus & arrayHasValues) )
		||   ( ! (arrayStatus & arrayHasDimensions) )
		||   ( value_ptr == (cfiScalar *)NULL ) )
		{
			FAILURE_MSG("cannot get array values");
			return (cfiFailure);
		}
	}

	/*
	 * create a matrix from which the values in the subarray
	 * can be retrieved, via the function getMatrixValues()
	 */

	axisCount = getArrayAxisCount(array);

	colCount = 1;
	for ( i = 0; i < axisCount; i++ )
	{
		if ( indices[i] >= 0 )
		{
			subStart[i] = indices[i];
			colCount *= getArrayDimension(array, i);
		}
/* LINTED */
		else if ( indices[i] == UNDEFINED_COUNT )
		{
			break;
		}
		else
		{
			FAILURE_MSG("cannot get array values");
			return (cfiFailure);
		}
	}

	rowCount = 1;
	for ( ; i < axisCount; i++ )
	{
		/*
		 * check for one (and only one) contiguous set
		 * of UNDEFINED_COUNT array index values
		 */

/* LINTED */
		if ( indices[i] == UNDEFINED_COUNT )
		{
			subStart[i] = 0;
			rowCount *= getArrayDimension(array, i);
		}
		else if ( indices[i] >= 0 )
		{
			break;
		}
		else
		{
			FAILURE_MSG("cannot get array values");
			return (cfiFailure);
		}
	}

	for ( ; i < axisCount; i++ )
	{
		if ( indices[i] >= 0 )
		{
			subStart[i] = indices[i];
		}
		else
		{
			FAILURE_MSG("cannot get array values");
			return (cfiFailure);
		}
	}

	if ( (subMatrix = allocateMatrix(getArrayValueType(array), 0, 0))
		== (cfiMatrix *)NULL )
	{
		FAILURE_MSG("cannot get array values");
		return (cfiFailure);
	}

	errorStatus = cfiSuccess;

	errorStatus |= setMatrixDimension((cfiDataAxis)rowAxis, rowCount, subMatrix);
	errorStatus |= setMatrixDimension((cfiDataAxis)columnAxis, colCount, subMatrix);
	errorStatus |= setMatrixValuePtr(getArrayValuePtr(array, subStart), subMatrix);

	/*
	 * get the subarray values, which are in the first column
	 * of the matrix
	 */

/* LINTED */
	errorStatus |= getMatrixValues(subMatrix, UNDEFINED_COUNT, 0, value_ptr);

	/*
	 * clean up
	 */

	errorStatus |= setMatrixValuePtr((void *)NULL, subMatrix);
	deallocateMatrix(&subMatrix);

	if ( (errorStatus & cfiFailure) )
	{
		FAILURE_MSG("cannot get array values");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of getArrayValues() */



/************************************************************************
 *									*
 * setArrayValues							*
 * --------------							*
 *	Converts a contiguous block of scalars to the type of		*
 *	an array, and stores the resulting values in the array.		*
 *									*
 * NOTE-Passing a contiguous set of UNDEFINED_COUNT values		*
 *	in "indices" sets all of the values in a subarray.		*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
setArrayValues
(
	const cfiScalar  *value_ptr,		/* values to be stored       */
	cfiArray         *const array,		/* the destination array     */
	const cfiCounter *const indices		/* array indices of interest */
)

{
	cfiErrorStatus errorStatus;	/* the return value                  */

	cfiArrayStatus arrayStatus;	/* status of the array               */

	cfiCounter axisCount;		/* number of array axes              */
	cfiCounter subStart[MAX_AXIS];	/* indices for start of subarray     */

	cfiMatrix  *subMatrix;		/* simplifies access to subarray     */
	cfiCounter  rowCount;		/* number of rows in submatrix       */
	cfiCounter  colCount;		/* number of columns in submatrix    */
	cfiCounter  i;			/* loop control                      */

	if ( (g_debugStatus & cfiArgDebug) )
	{
		arrayStatus = checkArray(array);

		if ( ( ! (arrayStatus & arrayHasValues) )
		||   ( ! (arrayStatus & arrayHasDimensions) )
		||   ( value_ptr == (cfiScalar *)NULL ) )
		{
			FAILURE_MSG("cannot set array values");
			return (cfiFailure);
		}
	}

	/*
	 * create a matrix in which the values in the subarray
	 * can be set, via the function setMatrixValues()
	 */

	axisCount = getArrayAxisCount(array);

	colCount = 1;
	for ( i = 0; i < axisCount; i++ )
	{
		if ( indices[i] >= 0 )
		{
			subStart[i] = indices[i];
			colCount *= getArrayDimension(array, i);
		}
/* LINTED */
		else if ( indices[i] == UNDEFINED_COUNT )
		{
			break;
		}
		else
		{
			FAILURE_MSG("cannot set array values");
			return (cfiFailure);
		}
	}

	rowCount = 1;
	for ( ; i < axisCount; i++ )
	{
		/*
		 * check for one (and only one) contiguous set
		 * of UNDEFINED_COUNT array index values
		 */

/* LINTED */
		if ( indices[i] == UNDEFINED_COUNT )
		{
			subStart[i] = 0;
			rowCount *= getArrayDimension(array, i);
		}
		else if ( indices[i] >= 0 )
		{
			break;
		}
		else
		{
			FAILURE_MSG("cannot set array values");
			return (cfiFailure);
		}
	}

	for ( ; i < axisCount; i++ )
	{
		if ( indices[i] >= 0 )
		{
			subStart[i] = indices[i];
		}
		else
		{
			FAILURE_MSG("cannot set array values");
			return (cfiFailure);
		}
	}

	if ( (subMatrix = allocateMatrix(getArrayValueType(array), 0, 0))
		== (cfiMatrix *)NULL )
	{
		FAILURE_MSG("cannot set array values");
		return (cfiFailure);
	}

	errorStatus = cfiSuccess;

	errorStatus |= setMatrixDimension((cfiDataAxis)rowAxis, rowCount, subMatrix);
	errorStatus |= setMatrixDimension((cfiDataAxis)columnAxis, colCount, subMatrix);
	errorStatus |= setMatrixValuePtr(getArrayValuePtr(array, subStart), subMatrix);

	/*
	 * set the subarray values, which are in the first column
	 * of the matrix
	 */

/* LINTED */
	errorStatus |= setMatrixValues(value_ptr, subMatrix, UNDEFINED_COUNT, 0);

	/*
	 * clean up
	 */

	errorStatus |= setMatrixValuePtr((void *)NULL, subMatrix);
	deallocateMatrix(&subMatrix);

	if ( (errorStatus & cfiFailure) )
	{
		FAILURE_MSG("cannot set array values");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of setArrayValues() */



/************************************************************************
 *									*
 * allocateArray							*
 * -------------							*
 *	Returns a pointer to an array of the specified type		*
 *	and dimensions, with its values set to zero.			*
 *									*
 * NOTE-Passing a NULL "dimensions" pointer is equivalent		*
 *	to specifying that there are zero elements along		*
 *	each array axis.						*
 *									*
 ********************************************************************bwr*/

cfiArray *
allocateArray
(
	const cfiDataType  valueType,		/* type of array values     */
	const cfiCounter   axisCount,		/* number of array axes     */
	const cfiCounter  *const dimensions	/* #'s of values along axes */
)

{
	cfiArray *array;		/* the return value                 */

	cfiCounter dimProd;		/* product of the array dimensions  */
	cfiCounter dimension;		/* number of values along an axis   */
	cfiCounter i;			/* loop control                     */

	(void)sprintf(l_msg, "allocating a %ld-D array", axisCount);
	DEBUG_MSG(cfiMemoryDebug, l_msg);

	if ( ( axisCount < 0 )
	||   ( axisCount > MAX_AXIS )
/* LINTED */
	||   ( (array = (cfiArray *)malloc(sizeof(cfiArray)))
		== (cfiArray *)NULL ) )
	{
		FAILURE_MSG("cannot allocate array");
		return ((cfiArray *)NULL);
	}

	array->axisCount = axisCount;

	if ( axisCount == 0 )
	{
		dimProd = 0;
	}
	else
	{
		dimProd = 1;
		for ( i = 0; i < axisCount; i++ )
		{
			if ( dimensions != (cfiCounter *)NULL )
			{
				dimension = dimensions[i];
			}
			else
			{
				dimension = 0;
			}

			if ( i == 0 )
			{
				(void)sprintf(l_msg, "...with dimensions %ld", dimension);
			}
			else
			{
				(void)sprintf(l_msg, "%sx%ld", l_msg, dimension);
			}

			if ( dimension < 0 )
			{
				(void)sprintf(l_msg, "%s...", l_msg);
				DEBUG_MSG(cfiMemoryDebug, l_msg);

				FAILURE_MSG("cannot allocate array");
				deallocateArray(&array);
				return ((cfiArray*)NULL);
			}

			(array->dimensions)[i] = dimension;
			dimProd *= dimension;

		} /* end of for ( i = 0; ... ) loop */

		DEBUG_MSG(cfiMemoryDebug, l_msg);

	} /* end of if ( axisCount == 0 ) ... else ... statement */

	for ( i = axisCount; i < MAX_AXIS; i++ )
	{
		(array->dimensions)[i] = 0;
	}

	if ( (array->arrayValues = allocateVector(valueType, dimProd))
		== (cfiVector *)NULL )
	{
		FAILURE_MSG("cannot allocate array");
		deallocateArray(&array);
		return ((cfiArray *)NULL);
	}

	return (array);

} /* end of allocateArray() */



/************************************************************************
 *									*
 * deallocateArray							*
 * ---------------							*
 *	Deallocates memory for an array.				*
 *									*
 ********************************************************************bwr*/

void
deallocateArray
(
	cfiArray **const array_ptr	/* the array of interest */
)

{
	DEBUG_MSG(cfiMemoryDebug, "deallocating array");

	if ( (*array_ptr) == (cfiArray *)NULL )
	{
		return;
	}

	deallocateVector(&((*array_ptr)->arrayValues));
	(void)free(*array_ptr);
	(*array_ptr) = (cfiArray *)NULL;
	
} /* end of deallocateArray() */



/************************************************************************
 *									*
 * arrayOpArray								*
 * ------------								*
 *	Performs pointwise arithmetic operations on two arrays.		*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
arrayOpArray
(
	const cfiArray *const array1,		/* one array            */
	const cfiOp     op,			/* the operation        */
	const cfiArray *const array2,		/* the other array      */
	cfiArray       *const resultArray	/* holds the result     */
)

{
	cfiErrorStatus errorStatus;		/* the return value     */

	cfiArrayStatus arrayStatus;		/* status of the arrays */

	cfiCounter i;				/* loop control         */

	DEBUG_MSG(cfiTraceDebug, "operating pointwise on arrays");

	if ( (g_debugStatus & cfiArgDebug) )
	{
		arrayStatus = (checkArray(array1) & checkArray(array2) & checkArray(resultArray));

		if ( ( ! (arrayStatus & arrayHasValues) )
		||   ( ! (arrayStatus & arrayHasDimensions) )
		||   ( resultArray->axisCount != array1->axisCount )
		||   ( resultArray->axisCount != array2->axisCount ) )
		{
			FAILURE_MSG("cannot operate pointwise on arrays");
			return (cfiFailure);
		}

		for ( i = 0; i < resultArray->axisCount; i++ )
		{
			if ( ( (resultArray->dimensions)[i] != (array1->dimensions)[i] )
			||   ( (resultArray->dimensions)[i] != (array2->dimensions)[i] ) )
			{
				FAILURE_MSG("cannot operate pointwise on arrays");
				return (cfiFailure);
			}
		}
	}

	if ( ((errorStatus = vectorOpVector(array1->arrayValues, op, array2->arrayValues, resultArray->arrayValues))
		& cfiFailure) )
	{
		FAILURE_MSG("cannot operate pointwise on arrays");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of arrayOpArray() */



/************************************************************************
 *									*
 * scaleArray								*
 * ----------								*
 *	Multiplies an array by a scalar.				*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
scaleArray
(
	const cfiScalar  scalar,	/* the scalar            */
	cfiArray        *const array	/* the array of interest */
)

{
	cfiErrorStatus errorStatus;	/* the return value      */

	cfiArrayStatus arrayStatus;	/* status of the array   */

	DEBUG_MSG(cfiTraceDebug, "scaling array");

	if ( (g_debugStatus & cfiArgDebug) )
	{
		arrayStatus = checkArray(array);

		if ( ( ! (arrayStatus & arrayHasValues) )
		||   ( ! (arrayStatus & arrayHasDimensions) ) )
		{
			FAILURE_MSG("cannot scale array");
			return (cfiFailure);
		}
	}

	if ( ((errorStatus = scaleVector(scalar, array->arrayValues))
		& cfiFailure) )
	{
		FAILURE_MSG("cannot scale array");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of scaleArray() */



/************************************************************************
 *									*
 * copyArray								*
 * ----------								*
 *	Copies values from an array of arbitrary type			*
 *	to another array of arbitrary type.				*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
copyArray
(
	const cfiArray *const fromArray,	/* the source array      */
	cfiArray       *const toArray		/* the destination array */
)

{
	cfiErrorStatus errorStatus;		/* the return value      */

	cfiArrayStatus arrayStatus;		/* status of the arrays  */

	cfiCounter i;				/* loop control          */

	DEBUG_MSG(cfiTraceDebug, "copying array");

	if ( (g_debugStatus & cfiArgDebug) )
	{
		arrayStatus = (checkArray(fromArray) & checkArray(toArray));

		if ( ( ! (arrayStatus & arrayHasValues) )
		||   ( ! (arrayStatus & arrayHasDimensions) )
		||   ( fromArray->axisCount != toArray->axisCount ) )
		{
			FAILURE_MSG("cannot copy array");
			return (cfiFailure);
		}

		for ( i = 0; i < fromArray->axisCount; i++ )
		{
			if ( (fromArray->dimensions)[i] != (toArray->dimensions)[i] )
			{
				FAILURE_MSG("cannot copy array");
				return (cfiFailure);
			}
		}
	}

	if ( ((errorStatus = copyVector(fromArray->arrayValues, toArray->arrayValues))
		& cfiFailure) )
	{
		FAILURE_MSG("cannot copy array");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of copyArray() */



/************************************************************************
 *									*
 * writeArray								*
 * ----------								*
 *	Writes an array to a stream.					*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
writeArray
(
	const cfiArray *const array,	/* the array of interest */
	FILE           *const stream	/* the output stream     */
)

{
	cfiErrorStatus errorStatus;	/* the return value      */

	cfiCounter i;			/* loop control          */

	DEBUG_MSG(cfiInOutDebug, "writing array");

	if ( array == (cfiArray *)NULL )
	{
		if ( fprintf(stream, "null array\n") < 0 )
		{
			FAILURE_MSG("cannot write array");
			return (cfiFailure);
		}

		return (cfiSuccess);
	}

	(void)fprintf(stream, "array dimensions [\n");
	for ( i = 0; i < array->axisCount; i++ )
	{
		(void)fprintf(stream, " %ld", (array->dimensions)[i]);
	}
	(void) fprintf(stream, "\n]\n");

	if ( ((errorStatus = writeVector(array->arrayValues, (array->dimensions)[0], stream))
		& cfiFailure) )
	{
		FAILURE_MSG("cannot write array");
		return (cfiFailure);
	}
		
	return (errorStatus);

} /* end of writeArray() */



/************************************************************************
 *									*
 * checkArray								*
 * ----------								*
 *	Checks if an array has values of a known type			*
 *	and valid dimensions.						*
 *									*
 ********************************************************************bwr*/

cfiArrayStatus
checkArray
(
	const cfiArray *const array	/* the array of interest           */
)

{
	cfiArrayStatus arrayStatus;	/* the return value                */

	cfiCounter dimProd;		/* product of the array dimensions */
	cfiCounter i;			/* loop control                    */

	DEBUG_MSG(cfiTraceDebug, "checking array");

	arrayStatus = emptyArray;

	if ( array == (cfiArray *)NULL )
	{
		return (arrayStatus);
	}

	/*
	 * check for values of a known type
	 */

	if ( (checkVector(array->arrayValues) & vectorHasValues) )
	{
		arrayStatus |= arrayHasValues;
	}

	/*
	 * check for valid dimensions
	 */

	if ( array->axisCount == 0 )
	{
		return (arrayStatus);
	}

	dimProd = 1;
	for ( i = 0; i < array->axisCount; i++ )
	{
		if ( (dimProd *= (array->dimensions)[i])
			<= 0 )
		{
			return (arrayStatus);
		}
	}

	if ( dimProd == getVectorValueCount(array->arrayValues) )
	{
		arrayStatus |= arrayHasDimensions;
	}

	return (arrayStatus);

} /* end of checkArray() */



/************************************************************************
 *									*
 * zeroArray								*
 * ---------								*
 *	Sets the values in an array to zero.				*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
zeroArray
(
	cfiArray *const array		/* the array of interest   */
)

{
	cfiErrorStatus errorStatus;	/* the return value        */

	cfiArrayStatus arrayStatus;	/* the status of the array */

	DEBUG_MSG(cfiTraceDebug, "setting array to zero");

	if ( (g_debugStatus & cfiArgDebug) )
	{
		arrayStatus = checkArray(array);

		if ( ( ! (arrayStatus & arrayHasValues) )
		||   ( ! (arrayStatus & arrayHasDimensions) ) )
		{
			FAILURE_MSG("cannot set array to zero");
			return (cfiFailure);
		}
	}

	if ( ((errorStatus = zeroVector(array->arrayValues))
		& cfiFailure) )
	{
		FAILURE_MSG("cannot set array to zero");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of zeroArray() */



/************************************************************************
 *									*
 * getArrayValueType							*
 * -----------------							*
 *	Returns the type of values in an array.				*
 *									*
 ********************************************************************bwr*/

cfiDataType
getArrayValueType
(
	const cfiArray *const array	/* the array of interest */
)

{
	cfiDataType valueType;		/* the return value      */

	if ( ( array == (cfiArray *)NULL )
	||   ( (valueType = getVectorValueType(array->arrayValues))
		== unknownType ) )
	{
		DEBUG_MSG(cfiReturnDebug, "type of array values is unknown");
		return (unknownType);
	}

	return (valueType);

} /* end of getArrayValueType() */



/************************************************************************
 *									*
 * setArrayValueType							*
 * -----------------							*
 *	Sets the type of values in an array.				*
 *									*
 * NOTE-To maintain a minimum level of data integrity,			*
 *	the array must have a NULL value pointer.			*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
setArrayValueType
(
	const cfiDataType  valueType,	/* the desired type of values */
	cfiArray          *const array	/* the array of interest      */
)

{
	cfiErrorStatus errorStatus;	/* the return value           */

	if ( ( array == (cfiArray *)NULL )
	||   ( ((errorStatus = setVectorValueType(valueType, array->arrayValues))
		& cfiFailure) ) )
	{
		FAILURE_MSG("cannot set type of array values");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of setArrayValueType() */



/************************************************************************
 *									*
 * getArrayAxisCount							*
 * -----------------							*
 *	Returns the number of axes for an array.			*
 *									*
 ********************************************************************bwr*/

cfiCounter
getArrayAxisCount
(
	const cfiArray *const array	/* the array of interest */
)

{
	if ( ( array == (cfiArray *)NULL )
	||   ( array->axisCount < 0 )
	||   ( array->axisCount > MAX_AXIS ) )
	{
		DEBUG_MSG(cfiReturnDebug, "number of array axes is ill-defined");
		return (0);
	}

	return (array->axisCount);

} /* end of getArrayAxisCount() */



/************************************************************************
 *									*
 * setArrayAxisCount							*
 * -----------------							*
 *	Sets the number of axes for an array.				*
 *									*
 * NOTE-To maintain a minimum level of data integrity,			*
 *	the array must have a NULL value pointer.			*
 *									*
 *	The dimensions for added or deleted axes			*
 *	are set to zero.						*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
setArrayAxisCount
(
	const cfiCounter  axisCount,	/* the desired number of axes      */
	cfiArray         *const array	/* the array of interest           */
)

{
	cfiErrorStatus errorStatus;	/* the return value                */

	cfiCounter dimProd;		/* product of the array dimensions */
	cfiCounter i;			/* loop control                    */

	if ( ( array == (cfiArray *)NULL )
	||   ( axisCount < 0 )
	||   ( axisCount > MAX_AXIS ) )
	{
		FAILURE_MSG("cannot set number of array axes");
		return (cfiFailure);
	}

	dimProd = 1;
	for ( i = 0; i < axisCount; i++ )
	{
		if ( i >= array->axisCount )
		{
			dimProd = 0;		/* (adding axes) */
			break;
		}

		dimProd *= (array->dimensions)[i];
	}

	/*
	 * set the dimensions for the added or deleted axes
	 */

	for ( i = array->axisCount; i < axisCount; i++ )
	{
		(array->dimensions)[i] = 0;	/* (added axis) */
	}

	for ( i = axisCount; i < array->axisCount; i++ )
	{
		(array->dimensions)[i] = 0;	/* (deleted axis) */
	}

	/*
	 * set the total storage count to be consistent
	 * with the array dimensions
	 */

	array->axisCount = axisCount;

	if ( ((errorStatus = setVectorValueCount(dimProd, array->arrayValues))
		& cfiFailure) )
	{
		FAILURE_MSG("cannot set number of array axes");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of setArrayAxisCount() */



/************************************************************************
 *									*
 * getArrayDimension							*
 * -----------------							*
 *	Returns the number of values along an array axis.		*
 *									*
 ********************************************************************bwr*/

cfiCounter
getArrayDimension
(
	const cfiArray    *const array,		/* the array of interest */
	const cfiDataAxis  axis			/* the axis of interest  */
)

{
	if ( ( array == (cfiArray *)NULL )
	||   ( axis < 0 )
	||   ( axis >= array->axisCount )
	||   ( (array->dimensions)[axis] <= 0 ) )
	{
		DEBUG_MSG(cfiReturnDebug, "array dimension is ill-defined");
		return (0);
	}

	return ((array->dimensions)[axis]);

} /* end of getArrayDimension() */



/************************************************************************
 *									*
 * setArrayDimension							*
 * -----------------							*
 *	Sets the number of values along an array axis.			*
 *									*
 * NOTE-To maintain a minimum level of data integrity,			*
 *	the array must have a NULL value pointer.			*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
setArrayDimension
(
	const cfiDataAxis  axis,	/* the axis of interest            */
	const cfiCounter   dimension,	/* number of values along axis     */
	cfiArray          *const array	/* the array of interest           */
)

{
	cfiErrorStatus errorStatus;	/* the return value                */

	cfiCounter dimProd;		/* product of the array dimensions */
	cfiCounter i;			/* loop control                    */

	if ( ( array == (cfiArray *)NULL )
	||   ( axis < 0 )
	||   ( axis >= array->axisCount )
	||   ( dimension < 0 ) )
	{
		FAILURE_MSG("cannot set array dimension");
		return (cfiFailure);
	}

	(array->dimensions)[axis] = dimension;

	/*
	 * set the total storage count to be consistent
	 * with the array dimensions
	 */

	dimProd = 1;
	for ( i = 0; i < array->axisCount; i++ )
	{
		dimProd *= (array->dimensions)[i];
	}

	if ( ((errorStatus = setVectorValueCount(dimProd, array->arrayValues))
		& cfiFailure) )
	{
		FAILURE_MSG("cannot set array dimension");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of setArrayDimension() */



/************************************************************************
 *									*
 * getArrayValuePtr							*
 * ----------------							*
 *	Returns a void pointer to a value in an array.			*
 *									*
 * NOTE-Passing a NULL "indices" pointer is equivalent			*
 *	to specifying the [0][0]...[0]-th element			*
 *	of the array as the value of interest.				*
 *									*
 ********************************************************************bwr*/

void *
getArrayValuePtr
(
	const cfiArray   *const array,	/* the array of interest              */
	const cfiCounter *const indices	/* array indices of value of interest */
)

{
	void *value_ptr;		/* the return value                   */

	cfiCounter index;		/* vector index for an array value    */
	cfiCounter stride;		/* index increment between values     */
	cfiCounter i;			/* loop control                       */

	if ( array == (cfiArray *)NULL )
	{
		DEBUG_MSG(cfiReturnDebug, "cannot get pointer to array value");
		return ((void *)NULL);
	}
	
	index = 0;	/* start at the [0][0]...[0]-th array element */

	if ( indices != (cfiCounter *)NULL )
	{
		/*
		 * stride across columns, rows, planes, etc.
		 * of array to reach the value of interest
		 */

		stride = 1;
		for ( i = 0; i < array->axisCount; i++ )
		{
/* LINTED */
			if ( indices[i] == UNDEFINED_COUNT )
			{
				DEBUG_MSG(cfiReturnDebug, "cannot get pointer to array value");
				return ((void *)NULL);
			}

			index += indices[i] * stride;
			stride *= (array->dimensions)[i];
		}
	}

	if ( (value_ptr = getVectorValuePtr(array->arrayValues, index))
		== (void *)NULL )
	{
		DEBUG_MSG(cfiReturnDebug, "cannot get pointer to array value");
		return ((void *)NULL);
	}

	return (value_ptr);

} /* end of getArrayValuePtr() */



/************************************************************************
 *									*
 * setArrayValuePtr							*
 * ----------------							*
 *	Sets the pointer to the values in an array.			*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
setArrayValuePtr
(
	void     *const value_ptr,	/* pointer to values of interest */
	cfiArray *const array		/* the array of interest         */
)

{
	cfiErrorStatus errorStatus;	/* the return value              */

	if ( ( array == (cfiArray *)NULL )
	||   ( ((errorStatus = setVectorValuePtr(value_ptr, array->arrayValues))
		& cfiFailure) ) )
	{
		FAILURE_MSG("cannot set pointer to array values");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of setArrayValuePtr() */
