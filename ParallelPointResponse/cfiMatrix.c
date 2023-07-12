/*
 * SCCS:  @(#)cfiMatrix.c  1.5  6/7/95  11:33:48
 */

/************************************************************************
 *									*
 * -----------------							*
 * FUNCTIONS HEREIN:							*
 * -----------------							*
 *									*
 * matrixTimesMatrix							*
 * -----------------							*
 *	Multiplies two matrices of arbitrary types and			*
 *	stores the result in a matrix of arbitrary type.		*
 *									*
 * matrixTimesVector							*
 * -----------------							*
 *	Multiplies a matrix and a column vector of arbitrary types	*
 *	and stores the result in a vector of arbitrary type.		*
 *									*
 * vectorTimesMatrix							*
 * -----------------							*
 *	Multiplies a row vector and a matrix of arbitrary types		*
 *	and stores the result in a vector of arbitrary type.		*
 *									*
 * allocateMatrix							*
 * --------------							*
 *	Returns a pointer to a matrix of the specified type and		*
 *	numbers of rows and columns, with its values set to zero.	*
 *									*
 * deallocateMatrix							*
 * ----------------							*
 *	Deallocates memory for a matrix.				*
 *									*
 * matrixOpMatrix							*
 * --------------							*
 *	Performs pointwise arithmetic operations on two matrices.	*
 *									*
 * scaleMatrix								*
 * -----------								*
 *	Multiplies a matrix by a scalar.				*
 *									*
 * copyMatrix								*
 * ----------								*
 *	Copies values from a matrix of arbitrary type			*
 *	to another matrix of arbitrary type.				*
 *									*
 * writeMatrix								*
 * -----------								*
 *	Writes a matrix to a stream.					*
 *									*
 * checkMatrix								*
 * -----------								*
 *	Checks if a matrix has values of a known type			*
 *	and valid dimensions.						*
 *									*
 * zeroMatrix								*
 * ----------								*
 *	Sets the values in a matrix to zero.				*
 *									*
 * getMatrixValues							*
 * ---------------							*
 *	Retrieves values from a matrix, and stores the values		*
 *	in a contiguous block of scalars.				*
 *									*
 * setMatrixValues							*
 * ---------------							*
 *	Converts a contiguous block of scalars to the type of		*
 *	a matrix, and stores the resulting values in the matrix.	*
 *									*
 * getMatrixValueType							*
 * ------------------							*
 *	Returns the type of values in a matrix.				*
 *									*
 * setMatrixValueType							*
 * ------------------							*
 *	Sets the type of values in a matrix.				*
 *									*
 * getMatrixDimension							*
 * ------------------							*
 *	Returns the number of rows or columns in a matrix.		*
 *									*
 * setMatrixDimension							*
 * ------------------							*
 *	Sets the number of rows or columns in a matrix.			*
 *									*
 * getMatrixValuePtr							*
 * -----------------							*
 *	Returns a void pointer to a value in a matrix.			*
 *									*
 * setMatrixValuePtr							*
 * -----------------							*
 *	Sets the pointer to the values in a matrix.			*
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



static char l_msg[MSG_LENGTH];		/* for compound messages */



/************************************************************************
 *									*
 * matrixTimesMatrix							*
 * -----------------							*
 *	Multiplies two matrices of arbitrary types and			*
 *	stores the result in a matrix of arbitrary type.		*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
matrixTimesMatrix
(
	const cfiMatrix *const matrix1,		/* one matrix             */
	const cfiMatrix *const matrix2,		/* the other matrix       */
	cfiMatrix       *const productMatrix	/* their product          */
)

{
	cfiErrorStatus errorStatus;		/* the return value       */

	cfiMatrixStatus matrixStatus;		/* status of the matrices */

	static cfiCounter  vecBufferLength = 0;	/* for internal buffers   */
	static cfiCounter  outBufferLength = 0;
	static cfiVector  *rowBuffer = (cfiVector *)NULL;
	static cfiVector  *colBuffer = (cfiVector *)NULL;
	static cfiMatrix  *outBuffer = (cfiMatrix *)NULL;
	static cfiScalar  *rowValues = (cfiScalar *)NULL;
	static cfiScalar  *colValues = (cfiScalar *)NULL;
	static cfiScalar  *outValues = (cfiScalar *)NULL;
	cfiScalar         *value_ptr;

	cfiCounter rowCount;			/* loop control           */
	cfiCounter colCount;
	cfiCounter innerCount;
	cfiCounter row; 
	cfiCounter col; 

	DEBUG_MSG(cfiTraceDebug, "multiplying matrices");

	rowCount = getMatrixDimension(matrix1, (cfiDataAxis)rowAxis);
	colCount = getMatrixDimension(matrix2, (cfiDataAxis)columnAxis);
	innerCount = getMatrixDimension(matrix1, (cfiDataAxis)columnAxis);

	if ( (g_debugStatus & cfiArgDebug) )
	{
		matrixStatus = (checkMatrix(matrix1) & checkMatrix(matrix2) & checkMatrix(productMatrix));

		if ( ( ! (matrixStatus & matrixHasValues) )
		||   ( ! (matrixStatus & matrixHasDimensions) )
		||   ( rowCount != getMatrixDimension(productMatrix, (cfiDataAxis)rowAxis) )
		||   ( colCount != getMatrixDimension(productMatrix, (cfiDataAxis)columnAxis) )
		||   ( innerCount != getMatrixDimension(matrix2, (cfiDataAxis)rowAxis) ) )
		{
			FAILURE_MSG("cannot multiply matrices");
			return (cfiFailure);
		}
	}

	/*
	 * set up the internal buffers
	 */

	if ( innerCount > vecBufferLength )
	{
		deallocateVector(&rowBuffer);
		deallocateVector(&colBuffer);

		if ( ( (rowBuffer = allocateVector((cfiDataType)scalarType, innerCount))
			== (cfiVector *)NULL )
		||   ( (colBuffer = allocateVector((cfiDataType)scalarType, innerCount))
			== (cfiVector *)NULL ) )
		{
			FAILURE_MSG("cannot multiply matrices");
			return (cfiFailure);
		}

		rowValues = (cfiScalar *)getVectorValuePtr(rowBuffer, 0);
		colValues = (cfiScalar *)getVectorValuePtr(colBuffer, 0);
		vecBufferLength = innerCount;
	}

	if ( (rowCount * colCount) > outBufferLength )
	{
		deallocateMatrix(&outBuffer);

		if ( (outBuffer = allocateMatrix((cfiDataType)scalarType, rowCount, colCount))
			== (cfiMatrix *)NULL )
		{
			FAILURE_MSG("cannot multiply matrices");
			return (cfiFailure);
		}

		outValues = (cfiScalar *)getMatrixValuePtr(outBuffer, 0, 0);
		outBufferLength = rowCount * colCount;
	}

	errorStatus = cfiSuccess;

	errorStatus |= setVectorValuePtr((void *)NULL, rowBuffer);
	errorStatus |= setVectorValueCount(innerCount, rowBuffer);
	errorStatus |= setVectorValuePtr((void *)rowValues, rowBuffer);

	errorStatus |= setVectorValuePtr((void *)NULL, colBuffer);
	errorStatus |= setVectorValueCount(innerCount, colBuffer);
	errorStatus |= setVectorValuePtr((void *)colValues, colBuffer);

	errorStatus |= setMatrixValuePtr((void *)NULL, outBuffer);
	errorStatus |= setMatrixDimension((cfiDataAxis)rowAxis, rowCount, outBuffer);
	errorStatus |= setMatrixDimension((cfiDataAxis)columnAxis, colCount, outBuffer);
	errorStatus |= setMatrixValuePtr((void *)outValues, outBuffer);

	/*
	 * multiply the matrices
	 */

	for ( row = 0; row < rowCount; row++ )
	{
/* LINTED */
		errorStatus |= getMatrixValues(matrix1, row, UNDEFINED_COUNT, rowValues);

		value_ptr = (cfiScalar *)getMatrixValuePtr(outBuffer, row, 0);

		for ( col = 0; col < colCount; col++ )
		{
/* LINTED */
			errorStatus |= getMatrixValues(matrix2, UNDEFINED_COUNT, col, colValues);
			errorStatus |= vectorDotVector(rowBuffer, colBuffer, value_ptr);

			value_ptr++;
		}
	}

	errorStatus |= copyMatrix(outBuffer, productMatrix);

	if ( (errorStatus & cfiFailure) )
	{
		FAILURE_MSG("cannot multiply matrices");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of matrixTimesMatrix() */



/************************************************************************
 *									*
 * matrixTimesVector							*
 * -----------------							*
 *	Multiplies a matrix and a column vector of arbitrary types	*
 *	and stores the result in a vector of arbitrary type.		*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
matrixTimesVector
(
	const cfiMatrix *const matrix,		/* the matrix            */
	const cfiVector *const colVector,	/* the column vector     */
	cfiVector       *const productVector	/* their product         */
)

{
	cfiErrorStatus errorStatus;		/* the return value      */

	cfiMatrixStatus matrixStatus;		/* status of the matrix  */
	cfiVectorStatus vectorStatus;		/* status of the vectors */

	static cfiCounter  rowBufferLength = 0;	/* for internal buffers   */
	static cfiCounter  outBufferLength = 0;
	static cfiVector  *rowBuffer = (cfiVector *)NULL;
	static cfiVector  *outBuffer = (cfiVector *)NULL;
	static cfiScalar  *rowValues = (cfiScalar *)NULL;
	static cfiScalar  *outValues = (cfiScalar *)NULL;
	cfiScalar         *value_ptr;

	cfiCounter rowCount;			/* loop control          */
	cfiCounter innerCount;
	cfiCounter row; 

	DEBUG_MSG(cfiTraceDebug, "multiplying matrix and column vector");

	rowCount = getMatrixDimension(matrix, (cfiDataAxis)rowAxis);
	innerCount = getMatrixDimension(matrix, (cfiDataAxis)columnAxis);

	if ( (g_debugStatus & cfiArgDebug) )
	{
		matrixStatus = checkMatrix(matrix);
		vectorStatus = (checkVector(colVector) & checkVector(productVector));

		if ( ( ! (matrixStatus & matrixHasValues) )
		||   ( ! (matrixStatus & matrixHasDimensions) )
		||   ( ! (vectorStatus & vectorHasValues) )
		||   ( rowCount != getVectorValueCount(productVector) )
		||   ( innerCount != getVectorValueCount(colVector) ) )
		{
			FAILURE_MSG("cannot multiply matrix and column vector");
			return (cfiFailure);
		}
	}

	/*
	 * set up the internal buffers
	 */

	if ( innerCount > rowBufferLength )
	{
		deallocateVector(&rowBuffer);

		if ( (rowBuffer = allocateVector((cfiDataType)scalarType, innerCount))
			== (cfiVector *)NULL )
		{
			FAILURE_MSG("cannot multiply matrix and column vector");
			return (cfiFailure);
		}

		rowValues = (cfiScalar *)getVectorValuePtr(rowBuffer, 0);
		rowBufferLength = innerCount;
	}

	if ( rowCount > outBufferLength )
	{
		deallocateVector(&outBuffer);

		if ( (outBuffer = allocateVector((cfiDataType)scalarType, rowCount))
			== (cfiVector *)NULL )
		{
			FAILURE_MSG("cannot multiply matrix and column vector");
			return (cfiFailure);
		}

		outValues = (cfiScalar *)getVectorValuePtr(outBuffer, 0);
		outBufferLength = rowCount;
	}

	errorStatus = cfiSuccess;

	errorStatus |= setVectorValuePtr((void *)NULL, rowBuffer);
	errorStatus |= setVectorValueCount(innerCount, rowBuffer);
	errorStatus |= setVectorValuePtr((void *)rowValues, rowBuffer);

	errorStatus |= setVectorValuePtr((void *)NULL, outBuffer);
	errorStatus |= setVectorValueCount(rowCount, outBuffer);
	errorStatus |= setVectorValuePtr((void *)outValues, outBuffer);

	/*
	 * multiply the matrix and the column vector
	 */

	value_ptr = outValues;

	for ( row = 0; row < rowCount; row++ )
	{
/* LINTED */
		errorStatus |= getMatrixValues(matrix, row, UNDEFINED_COUNT, rowValues);
		errorStatus |= vectorDotVector(rowBuffer, colVector, value_ptr);

		value_ptr++;
	}

	errorStatus |= copyVector(outBuffer, productVector);

	if ( (errorStatus & cfiFailure) )
	{
		FAILURE_MSG("cannot multiply matrix and column vector");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of matrixTimesVector() */



/************************************************************************
 *									*
 * vectorTimesMatrix							*
 * -----------------							*
 *	Multiplies a row vector and a matrix of arbitrary types		*
 *	and stores the result in a vector of arbitrary type.		*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
vectorTimesMatrix
(
	const cfiVector *const rowVector,	/* the row vector        */
	const cfiMatrix *const matrix,		/* the matrix            */
	cfiVector       *const productVector	/* their product         */
)

{
	cfiErrorStatus errorStatus;		/* the return value      */

	cfiMatrixStatus matrixStatus;		/* status of the matrix  */
	cfiVectorStatus vectorStatus;		/* status of the vectors */

	static cfiCounter  colBufferLength = 0;	/* for internal buffers   */
	static cfiCounter  outBufferLength = 0;
	static cfiVector  *colBuffer = (cfiVector *)NULL;
	static cfiVector  *outBuffer = (cfiVector *)NULL;
	static cfiScalar  *colValues = (cfiScalar *)NULL;
	static cfiScalar  *outValues = (cfiScalar *)NULL;
	cfiScalar         *value_ptr;

	cfiCounter colCount;			/* loop control          */
	cfiCounter innerCount;
	cfiCounter col; 

	DEBUG_MSG(cfiTraceDebug, "multiplying row vector and matrix");

	colCount = getMatrixDimension(matrix, (cfiDataAxis)columnAxis);
	innerCount = getMatrixDimension(matrix, (cfiDataAxis)rowAxis);

	if ( (g_debugStatus & cfiArgDebug) )
	{
		matrixStatus = checkMatrix(matrix);
		vectorStatus = (checkVector(rowVector) & checkVector(productVector));

		if ( ( ! (matrixStatus & matrixHasValues) )
		||   ( ! (matrixStatus & matrixHasDimensions) )
		||   ( ! (vectorStatus & vectorHasValues) )
		||   ( colCount != getVectorValueCount(productVector) )
		||   ( innerCount != getVectorValueCount(rowVector) ) )
		{
			FAILURE_MSG("cannot multiply row vector and matrix");
			return (cfiFailure);
		}
	}

	/*
	 * set up the internal buffers
	 */

	if ( innerCount > colBufferLength )
	{
		deallocateVector(&colBuffer);

		if ( (colBuffer = allocateVector((cfiDataType)scalarType, innerCount))
			== (cfiVector *)NULL )
		{
			FAILURE_MSG("cannot multiply row vector and matrix");
			return (cfiFailure);
		}

		colValues = (cfiScalar *)getVectorValuePtr(colBuffer, 0);
		colBufferLength = innerCount;
	}

	if ( colCount > outBufferLength )
	{
		deallocateVector(&outBuffer);

		if ( (outBuffer = allocateVector((cfiDataType)scalarType, colCount))
			== (cfiVector *)NULL )
		{
			FAILURE_MSG("cannot multiply row vector and matrix");
			return (cfiFailure);
		}

		outValues = (cfiScalar *)getVectorValuePtr(outBuffer, 0);
		outBufferLength = colCount;
	}

	errorStatus = cfiSuccess;

	errorStatus |= setVectorValuePtr((void *)NULL, colBuffer);
	errorStatus |= setVectorValueCount(innerCount, colBuffer);
	errorStatus |= setVectorValuePtr((void *)colValues, colBuffer);

	errorStatus |= setVectorValuePtr((void *)NULL, outBuffer);
	errorStatus |= setVectorValueCount(colCount, outBuffer);
	errorStatus |= setVectorValuePtr((void *)outValues, outBuffer);

	/*
	 * multiply the row vector and the matrix
	 */

	value_ptr = outValues;

	for ( col = 0; col < colCount; col++ )
	{
/* LINTED */
		errorStatus |= getMatrixValues(matrix, UNDEFINED_COUNT, col, colValues);
		errorStatus |= vectorDotVector(rowVector, colBuffer, value_ptr);

		value_ptr++;
	}

	errorStatus |= copyVector(outBuffer, productVector);

	if ( (errorStatus & cfiFailure) )
	{
		FAILURE_MSG("cannot multiply row vector and matrix");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of vectorTimesMatrix() */



/************************************************************************
 *									*
 * allocateMatrix							*
 * --------------							*
 *	Returns a pointer to a matrix of the specified type and		*
 *	numbers of rows and columns, with its values set to zero.	*
 *									*
 ********************************************************************bwr*/

cfiMatrix *
allocateMatrix
(
	const cfiDataType valueType,	/* type of matrix values           */
	const cfiCounter  rowCount,	/* number of rows in the matrix    */
	const cfiCounter  colCount	/* number of columns in the matrix */
)

{
	cfiMatrix *matrix;		/* the return value                */

	(void)sprintf(l_msg, "allocating a %ldx%ld matrix", rowCount, colCount);
	DEBUG_MSG(cfiMemoryDebug, l_msg);

	if ( ( rowCount < 0 )
	||   ( colCount < 0 )
/* LINTED */
	||   ( (matrix = (cfiMatrix *)malloc(sizeof(cfiMatrix)))
		== (cfiMatrix *)NULL ) )
	{
		FAILURE_MSG("cannot allocate matrix");
		return ((cfiMatrix *)NULL);
	}
	
	matrix->rowCount = rowCount;
	matrix->colCount = colCount;

	if ( (matrix->matrixValues = allocateVector(valueType, rowCount * colCount))
		== (cfiVector *)NULL )
	{
		FAILURE_MSG("cannot allocate matrix");
		deallocateMatrix(&matrix);
		return ((cfiMatrix *)NULL);
	}

	return (matrix);

} /* end of allocateMatrix() */



/************************************************************************
 *									*
 * deallocateMatrix							*
 * ----------------							*
 *	Deallocates memory for a matrix.				*
 *									*
 ********************************************************************bwr*/

void
deallocateMatrix
(
	cfiMatrix **const matrix_ptr	/* the matrix of interest */
)

{
	DEBUG_MSG(cfiMemoryDebug, "deallocating matrix");

	if ( (*matrix_ptr) == (cfiMatrix *)NULL )
	{
		return;
	}

	deallocateVector(&((*matrix_ptr)->matrixValues));
	(void)free(*matrix_ptr);
	(*matrix_ptr) = (cfiMatrix *)NULL;
	
} /* end of deallocateMatrix() */



/************************************************************************
 *									*
 * matrixOpMatrix							*
 * --------------							*
 *	Performs pointwise arithmetic operations on two matrices.	*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
matrixOpMatrix
(
	const cfiMatrix *const matrix1,		/* one matrix             */
	const cfiOp      op,			/* the operation          */
	const cfiMatrix *const matrix2,		/* the other matrix       */
	cfiMatrix       *const resultMatrix	/* holds the result       */
)

{
	cfiErrorStatus errorStatus;		/* the return value       */

	cfiMatrixStatus matrixStatus;		/* status of the matrices */

	DEBUG_MSG(cfiTraceDebug, "operating pointwise on matrices");

	if ( (g_debugStatus & cfiArgDebug) )
	{
		matrixStatus = (checkMatrix(matrix1) & checkMatrix(matrix2) & checkMatrix(resultMatrix));

		if ( ( ! (matrixStatus & matrixHasValues) )
		||   ( ! (matrixStatus & matrixHasDimensions) )
		||   ( resultMatrix->rowCount != matrix1->rowCount )
		||   ( resultMatrix->rowCount != matrix2->rowCount )
		||   ( resultMatrix->colCount != matrix1->colCount )
		||   ( resultMatrix->colCount != matrix2->colCount ) )
		{
			FAILURE_MSG("cannot operate pointwise on matrices");
			return (cfiFailure);
		}
	}

	if ( ((errorStatus = vectorOpVector(matrix1->matrixValues, op, matrix2->matrixValues, resultMatrix->matrixValues))
		& cfiFailure) )
	{
		FAILURE_MSG("cannot operate pointwise on matrices");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of matrixOpMatrix() */



/************************************************************************
 *									*
 * scaleMatrix								*
 * -----------								*
 *	Multiplies a matrix by a scalar.				*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
scaleMatrix
(
	const cfiScalar  scalar,	/* the scalar             */
	cfiMatrix       *const matrix	/* the matrix of interest */
)

{
	cfiErrorStatus errorStatus;	/* the return value       */

	cfiMatrixStatus matrixStatus;	/* status of the matrix   */

	DEBUG_MSG(cfiTraceDebug, "scaling matrix");

	if ( (g_debugStatus & cfiArgDebug) )
	{
		matrixStatus = checkMatrix(matrix);

		if ( ( ! (matrixStatus & matrixHasValues) )
		||   ( ! (matrixStatus & matrixHasDimensions) ) )
		{
			FAILURE_MSG("cannot scale matrix");
			return (cfiFailure);
		}
	}

	if ( ((errorStatus = scaleVector(scalar, matrix->matrixValues))
		& cfiFailure) )
	{
		FAILURE_MSG("cannot scale matrix");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of scaleMatrix() */



/************************************************************************
 *									*
 * copyMatrix								*
 * ----------								*
 *	Copies values from a matrix of arbitrary type			*
 *	to another matrix of arbitrary type.				*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
copyMatrix
(
	const cfiMatrix *const fromMatrix,	/* the source matrix      */
	cfiMatrix       *const toMatrix		/* the destination matrix */
)

{
	cfiErrorStatus errorStatus;		/* the return value       */

	cfiMatrixStatus matrixStatus;		/* status of the matrices */

	DEBUG_MSG(cfiTraceDebug, "copying matrix");

	if ( (g_debugStatus & cfiArgDebug) )
	{
		matrixStatus = (checkMatrix(fromMatrix) & checkMatrix(toMatrix));

		if ( ( ! (matrixStatus & matrixHasValues) )
		||   ( ! (matrixStatus & matrixHasDimensions) )
		||   ( fromMatrix->rowCount != toMatrix->rowCount )
		||   ( fromMatrix->colCount != toMatrix->colCount ) )
		{
			FAILURE_MSG("cannot copy matrix");
			return (cfiFailure);
		}
	}

	if ( ((errorStatus = copyVector(fromMatrix->matrixValues, toMatrix->matrixValues))
		& cfiFailure) )
	{
		FAILURE_MSG("cannot copy matrix");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of copyMatrix() */



/************************************************************************
 *									*
 * writeMatrix								*
 * -----------								*
 *	Writes a matrix to a stream.					*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
writeMatrix
(
	const cfiMatrix *const matrix,	/* the matrix of interest */
	FILE            *const stream	/* the output stream      */
)

{
	cfiErrorStatus errorStatus;	/* the return value       */

	DEBUG_MSG(cfiInOutDebug, "writing matrix");

	if ( matrix == (cfiMatrix *)NULL )
	{
		if ( fprintf(stream, "null matrix\n") < 0 )
		{
			FAILURE_MSG("cannot write matrix");
			return (cfiFailure);
		}

		return (cfiSuccess);
	}

	(void)fprintf(stream, "matrix dimensions [\n");
	(void)fprintf(stream, " %ld %ld\n", matrix->rowCount, matrix->colCount);
	(void)fprintf(stream, "]\n");

	if ( ((errorStatus = writeVector(matrix->matrixValues, matrix->colCount, stream))
		& cfiFailure) )
	{
		FAILURE_MSG("cannot write matrix");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of writeMatrix() */



/************************************************************************
 *									*
 * checkMatrix								*
 * -----------								*
 *	Checks if a matrix has values of a known type			*
 *	and valid dimensions.						*
 *									*
 ********************************************************************bwr*/

cfiMatrixStatus
checkMatrix
(
	const cfiMatrix *const matrix	/* the matrix of interest */
)

{
	cfiMatrixStatus matrixStatus;	/* the return value       */

	DEBUG_MSG(cfiTraceDebug, "checking matrix");

	matrixStatus = emptyMatrix;

	if ( matrix == (cfiMatrix *)NULL )
	{
		return (matrixStatus);
	}

	/*
	 * check for values of a known type
	 */

	if ( (checkVector(matrix->matrixValues) & vectorHasValues) )
	{
		matrixStatus |= matrixHasValues;
	}

	/*
	 * check for valid dimensions
	 */

	if ( ( matrix->rowCount > 0 )
	&&   ( matrix->colCount > 0 )
	&&   ( (matrix->rowCount * matrix->colCount) == getVectorValueCount(matrix->matrixValues) ) )
	{
		matrixStatus |= matrixHasDimensions;
	}

	return (matrixStatus);

} /* end of checkMatrix() */



/************************************************************************
 *									*
 * zeroMatrix								*
 * ----------								*
 *	Sets the values in a matrix to zero.				*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
zeroMatrix
(
	cfiMatrix *const matrix		/* the matrix of interest   */
)

{
	cfiErrorStatus errorStatus;	/* the return value         */

	cfiMatrixStatus matrixStatus;	/* the status of the matrix */

	DEBUG_MSG(cfiTraceDebug, "setting matrix to zero"); 

	if ( (g_debugStatus & cfiArgDebug) )
	{
		matrixStatus = checkMatrix(matrix);

		if ( ( ! (matrixStatus & matrixHasValues) )
		||   ( ! (matrixStatus & matrixHasDimensions) ) )
		{
			FAILURE_MSG("cannot set matrix to zero");
			return (cfiFailure);
		}
	}

	if ( ((errorStatus = zeroVector(matrix->matrixValues))
		& cfiFailure) )
	{
		FAILURE_MSG("cannot set matrix to zero");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of zeroMatrix() */



/************************************************************************
 *									*
 * getMatrixValues							*
 * ---------------							*
 *	Retrieves values from a matrix, and stores the values		*
 *	in a contiguous block of scalars.				*
 *									*
 * NOTE-Passing a "row" index of UNDEFINED_COUNT gets all of		*
 *	the values from a column, while passing a "col" index		*
 *	of UNDEFINED_COUNT gets all of the values from a row.		*
 *	Setting both indices to UNDEFINED_COUNT gets all of		*
 *	the values from the matrix.					*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
getMatrixValues
(
	const cfiMatrix  *const matrix,	/* the source matrix                */
	const cfiCounter  row,		/* index of row of interest         */
	const cfiCounter  col,		/* index of column of interest      */
	cfiScalar        *value_ptr	/* destination of the values        */
)

{
	cfiMatrixStatus matrixStatus;	/* the status of the matrix         */

	void *matrixValue_ptr;		/* pointer to matrix values         */

	cfiCounter valueCount;		/* number of values to be retrieved */
	cfiCounter index;		/* vector index for a matrix value  */
	cfiCounter stride;		/* index increment between values   */

	cfiScalar     *end_value;	/* loop control                     */
	unsigned char *byte_ptr;
	short         *short_ptr;
	long          *long_ptr;
	float         *float_ptr;
	double        *double_ptr;

	matrixValue_ptr = getMatrixValuePtr(matrix, 0, 0);

	if ( (g_debugStatus & cfiArgDebug) )
	{
		matrixStatus = checkMatrix(matrix);

		if ( ( ! (matrixStatus & matrixHasValues) )
		||   ( ! (matrixStatus & matrixHasDimensions) )
		||   ( matrixValue_ptr == (void *)NULL )
		||   ( (row < 0) && (row != UNDEFINED_COUNT) )
		||   ( row >= matrix->rowCount )
		||   ( (col < 0) && (col != UNDEFINED_COUNT) )
		||   ( col >= matrix->colCount )
/* LINTED */
		||   ( value_ptr == (cfiScalar *)NULL ) )
		{
			FAILURE_MSG("cannot get matrix values");
			return (cfiFailure);
		}
	}

	index = 0;		/* start at the [0][0]-th matrix element */
	valueCount = 1;

	if ( row >= 0 )
	{
		/*
		 * getting values from just one row,
		 * so jump down to the start of the row
		 */

		index += row * matrix->colCount;
	}
	else
	{
		/*
		 * getting values from each row
		 */

		valueCount *= matrix->rowCount;
	}

	if ( col >= 0 )
	{
		/*
		 * getting values from just one column,
		 * so jump over to the column
		 */

		index += col;
		stride = matrix->colCount;
	}
	else
	{
		/*
		 * getting values from each column
		 */

		valueCount *= matrix->colCount;
		stride = 1;
	}

	end_value = value_ptr + valueCount;

	switch (getMatrixValueType(matrix))
	{
	  case byteType:

		byte_ptr = ((unsigned char *)matrixValue_ptr) + index;

		while ( value_ptr < end_value )
		{
			(*value_ptr++) = (cfiScalar)(*byte_ptr);
			byte_ptr += stride;
		}

		return (cfiSuccess);

	  case shortType:

		short_ptr = ((short *)matrixValue_ptr) + index;

		while ( value_ptr < end_value )
		{
			(*value_ptr++) = (cfiScalar)(*short_ptr);
			short_ptr += stride;
		}

		return (cfiSuccess);

	  case longType:

		long_ptr = ((long *)matrixValue_ptr) + index;

		while ( value_ptr < end_value )
		{
			(*value_ptr++) = (cfiScalar)(*long_ptr);
			long_ptr += stride;
		}

		return (cfiSuccess);

	  case floatType:

		float_ptr = ((float *)matrixValue_ptr) + index;

		while ( value_ptr < end_value )
		{
			(*value_ptr++) = (cfiScalar)(*float_ptr);
			float_ptr += stride;
		}

		return (cfiSuccess);

	  case doubleType:

		double_ptr = ((double *)matrixValue_ptr) + index;

		while ( value_ptr < end_value )
		{
			(*value_ptr++) = (cfiScalar)(*double_ptr);
			double_ptr += stride;
		}

		return (cfiSuccess);

	  default:

		break;

	} /* end of switch (getMatrixValueType(matrix)) */

	FAILURE_MSG("cannot get matrix values");
	return (cfiFailure);

} /* end of getMatrixValues() */



/************************************************************************
 *									*
 * setMatrixValues							*
 * ---------------							*
 *	Converts a contiguous block of scalars to the type of		*
 *	a matrix, and stores the resulting values in the matrix.	*
 *									*
 * NOTE-Passing a "row" index of UNDEFINED_COUNT sets all of		*
 *	the values in a column, while passing a "col" index		*
 *	of UNDEFINED_COUNT sets all of the values in a row.		*
 *	Setting both indices to UNDEFINED_COUNT sets all of		*
 *	the values in the matrix.					*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
setMatrixValues
(
	const cfiScalar  *value_ptr,	/* values to be stored in matrix   */
	cfiMatrix        *const matrix,	/* the destination matrix          */
	const cfiCounter  row,		/* index of row of interest        */
	const cfiCounter  col		/* index of column of interest     */
)

{
	cfiMatrixStatus matrixStatus;	/* the status of the matrix        */

	void *matrixValue_ptr;		/* pointer to matrix values        */

	cfiCounter valueCount;		/* number of values to be stored   */
	cfiCounter index;		/* linear index to a matrix value  */
	cfiCounter stride;		/* index increment between values  */

	const cfiScalar *end_value;	/* loop control                    */
	unsigned char   *byte_ptr;
	short           *short_ptr;
	long            *long_ptr;
	float           *float_ptr;
	double          *double_ptr;

	matrixValue_ptr = getMatrixValuePtr(matrix, 0, 0);

	if ( (g_debugStatus & cfiArgDebug) )
	{
		matrixStatus = checkMatrix(matrix);

		if ( ( ! (matrixStatus & matrixHasValues) )
		||   ( ! (matrixStatus & matrixHasDimensions) )
		||   ( matrixValue_ptr == (void *)NULL )
		||   ( (row < 0) && (row != UNDEFINED_COUNT) )
		||   ( row >= matrix->rowCount )
		||   ( (col < 0) && (col != UNDEFINED_COUNT) )
		||   ( col >= matrix->colCount )
/* LINTED */
		||   ( value_ptr == (cfiScalar *)NULL ) )
		{
			FAILURE_MSG("cannot set matrix values");
			return (cfiFailure);
		}
	}

	index = 0;		/* start at the [0][0]-th matrix element */
	valueCount = 1;

	if ( row >= 0 )
	{
		/*
		 * setting values in just one row,
		 * so jump down to the start of the row
		 */

		index += row * matrix->colCount;
	}
	else
	{
		/*
		 * setting values in each row
		 */

		valueCount *= matrix->rowCount;
	}

	if ( col >= 0 )
	{
		/*
		 * setting values in just one column,
		 * so jump over to the column
		 */

		index += col;
		stride = matrix->colCount;
	}
	else
	{
		/*
		 * setting values in each column
		 */

		valueCount *= matrix->colCount;
		stride = 1;
	}

	end_value = value_ptr + valueCount;

	switch (getMatrixValueType(matrix))
	{
	  case byteType:

		byte_ptr = ((unsigned char *)matrixValue_ptr) + index;

		while ( value_ptr < end_value )
		{
			if ( (*value_ptr) > 255.0 )
			{
				DEBUG_MSG(cfiMaxDebug, "clipping value of type byte");
				(*byte_ptr) = 255;
			}
			else if ( (*value_ptr) < 0.0 )
			{
				DEBUG_MSG(cfiMaxDebug, "clipping value of type byte");
				(*byte_ptr) = 0;
			}
			else
			{
				(*byte_ptr) = (unsigned char)((*value_ptr) + 0.5);
			}

			byte_ptr += stride;
			value_ptr++;

		} /* end of while ( value_ptr < end_value ) loop */

		return (cfiSuccess);

	  case shortType:

		short_ptr = ((short *)matrixValue_ptr) + index;

		while ( value_ptr < end_value )
		{
			if ( (*value_ptr) > (cfiScalar)MAXSHORT )
			{
				DEBUG_MSG(cfiMaxDebug, "clipping value of type short");
				(*short_ptr) = MAXSHORT;
			}
			else if ( (*value_ptr) < (cfiScalar)(-MAXSHORT) )
			{
				DEBUG_MSG(cfiMaxDebug, "clipping value of type short");
				(*short_ptr) = (-MAXSHORT);
			}
			else if ( (*value_ptr) >= 0.0 )
			{
				(*short_ptr) = (short)((*value_ptr) + 0.5);
			}
			else
			{
				(*short_ptr) = (short)((*value_ptr) - 0.5);
			}

			short_ptr += stride;
			value_ptr++;

		} /* end of while ( value_ptr < end_value ) loop */

		return (cfiSuccess);

	  case longType:

		long_ptr = ((long *)matrixValue_ptr) + index;

		while ( value_ptr < end_value )
		{
/* LINTED */
			if ( (*value_ptr) > (cfiScalar)MAXLONG )
			{
				DEBUG_MSG(cfiMaxDebug, "clipping value of type long");
/* LINTED */
				(*long_ptr) = MAXLONG;
			}
/* LINTED */
			else if ( (*value_ptr) < (cfiScalar)(-MAXLONG) )
			{
				DEBUG_MSG(cfiMaxDebug, "clipping value of type long");
/* LINTED */
				(*long_ptr) = (-MAXLONG);
			}
			else if ( (*value_ptr) >= 0.0 )
			{
				(*long_ptr) = (long)((*value_ptr) + 0.5);
			}
			else
			{
				(*long_ptr) = (long)((*value_ptr) - 0.5);
			}

			long_ptr += stride;
			value_ptr++;

		} /* end of while ( value_ptr < end_value ) loop */

		return (cfiSuccess);

	  case floatType:

		float_ptr = ((float *)matrixValue_ptr) + index;

		while ( value_ptr < end_value )
		{
			if ( (*value_ptr) > (cfiScalar)MAXFLOAT )
			{
				DEBUG_MSG(cfiMaxDebug, "clipping value of type float");
				(*float_ptr) = MAXFLOAT;
			}
			else if ( (*value_ptr) < (cfiScalar)(-MAXFLOAT) )
			{
				DEBUG_MSG(cfiMaxDebug, "clipping value of type float");
				(*float_ptr) = (-MAXFLOAT);
			}
			else
			{
				(*float_ptr) = (float)(*value_ptr);
			}

			float_ptr += stride;
			value_ptr++;

		} /* end of while ( value_ptr < end_value ) loop */

		return (cfiSuccess);

	  case doubleType:

		double_ptr = ((double *)matrixValue_ptr) + index;

		while ( value_ptr < end_value )
		{
			(*double_ptr) = (double)(*value_ptr);

			double_ptr += stride;
			value_ptr++;
	
		} /* end of while ( value_ptr < end_value ) loop */

		return (cfiSuccess);

	  default:

		break;

	} /* end of switch (getMatrixValueType(matrix)) */

	FAILURE_MSG("cannot set matrix values");
	return (cfiFailure);

} /* end of setMatrixValues() */



/************************************************************************
 *									*
 * getMatrixValueType							*
 * ------------------							*
 *	Returns the type of values in a matrix.				*
 *									*
 ********************************************************************bwr*/

cfiDataType
getMatrixValueType
(
	const cfiMatrix *const matrix	/* the matrix of interest */
)

{
	cfiDataType valueType;		/* the return value       */

	if ( ( matrix == (cfiMatrix *)NULL )
	||   ( (valueType = getVectorValueType(matrix->matrixValues))
		== unknownType ) )
	{
		DEBUG_MSG(cfiReturnDebug, "type of matrix values is unknown");
		return (unknownType);
	}

	return (valueType);

} /* end of getMatrixValueType() */



/************************************************************************
 *									*
 * setMatrixValueType							*
 * ------------------							*
 *	Sets the type of values in a matrix.				*
 *									*
 * NOTE-To maintain a minimum level of data integrity,			*
 *	the matrix must have a NULL value pointer.			*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
setMatrixValueType
(
	const cfiDataType  valueType,	/* the desired type of values */
	cfiMatrix         *const matrix	/* the matrix of interest     */
)

{
	cfiErrorStatus errorStatus;	/* the return value           */

	if ( ( matrix == (cfiMatrix *)NULL )
	||   ( ((errorStatus = setVectorValueType(valueType, matrix->matrixValues))
		& cfiFailure) ) )
	{
		FAILURE_MSG("cannot set type of matrix values");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of setMatrixValueType() */



/************************************************************************
 *									*
 * getMatrixDimension							*
 * ------------------							*
 *	Returns the number of rows or columns in a matrix.		*
 *									*
 ********************************************************************bwr*/

cfiCounter
getMatrixDimension
(
	const cfiMatrix   *const matrix,	/* the matrix of interest */
	const cfiDataAxis  axis			/* the axis of interest   */
)

{
	cfiCounter dimension;			/* the return value       */

	if ( matrix == (cfiMatrix *)NULL )
	{
		DEBUG_MSG(cfiReturnDebug, "matrix dimension is ill-defined");
		return (0);
	}

	if ( axis == (cfiDataAxis)columnAxis )
	{
		dimension = matrix->colCount;
	}
	else if ( axis == (cfiDataAxis)rowAxis )
	{
		dimension = matrix->rowCount;
	}
	else
	{
		DEBUG_MSG(cfiReturnDebug, "matrix dimension is ill-defined");
		return (0);
	}

	if ( dimension <= 0 )
	{
		DEBUG_MSG(cfiReturnDebug, "matrix dimension is ill-defined");
		return (0);
	}

	return (dimension);

} /* end of getMatrixDimension() */



/************************************************************************
 *									*
 * setMatrixDimension							*
 * ------------------							*
 *	Sets the number of rows or columns in a matrix.			*
 *									*
 * NOTE-To maintain a minimum level of data integrity,			*
 *	the matrix must have a NULL value pointer.			*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
setMatrixDimension
(
	const cfiDataAxis  axis,	/* the axis of interest   */
	const cfiCounter   dimension,	/* # of values along axis */
	cfiMatrix         *const matrix	/* the matrix of interest */
)

{
	cfiErrorStatus errorStatus;	/* the return value       */

	if ( ( matrix == (cfiMatrix *)NULL )
	||   ( dimension < 0 ) )
	{
		FAILURE_MSG("cannot set matrix dimension");
		return (cfiFailure);
	}

	if ( axis == (cfiDataAxis)columnAxis )
	{
		matrix->colCount = dimension;
	}
	else if ( axis == (cfiDataAxis)rowAxis )
	{
		matrix->rowCount = dimension;
	}
	else
	{
		FAILURE_MSG("cannot set matrix dimension");
		return (cfiFailure);
	}

	/*
	 * set the total storage count to be consistent
	 * with the matrix dimensions
	 */

	if ( ((errorStatus = setVectorValueCount(matrix->rowCount * matrix->colCount, matrix->matrixValues))
		& cfiFailure) )
	{
		FAILURE_MSG("cannot set matrix dimension");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of setMatrixDimension() */



/************************************************************************
 *									*
 * getMatrixValuePtr							*
 * -----------------							*
 *	Returns a void pointer to a value in a matrix.			*
 *									*
 ********************************************************************bwr*/

void *
getMatrixValuePtr
(
	const cfiMatrix  *const matrix,	/* the matrix of interest            */
	const cfiCounter  row,		/* row index of value of interest    */
	const cfiCounter  col		/* column index of value of interest */
)

{
	void *value_ptr;		/* the return value                  */

	if ( ( matrix == (cfiMatrix *)NULL )
	||   ( row == UNDEFINED_COUNT )
	||   ( col == UNDEFINED_COUNT )
	||   ( (value_ptr = getVectorValuePtr(matrix->matrixValues, (row * matrix->colCount) + col))
/* LINTED */
		== (void *)NULL ) )
	{
		DEBUG_MSG(cfiReturnDebug, "cannot get pointer to matrix value");
		return ((void *)NULL);
	}

	return (value_ptr);

} /* end of getMatrixValuePtr() */



/************************************************************************
 *									*
 * setMatrixValuePtr							*
 * -----------------							*
 *	Sets the pointer to the values in a matrix.			*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
setMatrixValuePtr
(
	void      *const value_ptr,	/* pointer to values of interest */
	cfiMatrix *const matrix		/* the matrix of interest        */
)

{
	cfiErrorStatus errorStatus;	/* the return value              */

	if ( ( matrix == (cfiMatrix *)NULL )
	||   ( ((errorStatus = setVectorValuePtr(value_ptr, matrix->matrixValues))
		& cfiFailure) ) )
	{
		FAILURE_MSG("cannot set pointer to matrix values");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of setMatrixValuePtr() */
