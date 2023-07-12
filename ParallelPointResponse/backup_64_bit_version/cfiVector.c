/*
 * SCCS:  @(#)cfiVector.c  1.9  6/7/95  11:33:50
 */

/************************************************************************
 *									*
 * -----------------							*
 * FUNCTIONS HEREIN:							*
 * -----------------							*
 *									*
 * vectorDotVector							*
 * ---------------							*
 *	Calculates the inner product of two vectors			*
 *	(accelerated for floating-point operands			*
 *	having the same type of values).				*
 *									*
 * vectorOpVector							*
 * --------------							*
 *	Performs pointwise arithmetic operations on two vectors		*
 *	(accelerated for floating-point operands and result		*
 *	having the same type of values).				*
 *									*
 * scaleVector								*
 * -----------								*
 *	Multiplies a vector by a scalar					*
 *	(accelerated for a floating-point vector).			*
 *									*
 * copyVector								*
 * ----------								*
 *	Copies values from a vector of arbitrary			*
 *	type to another vector of arbitrary type			*
 *	(accelerated for vectors in non-overlapping			*
 *	memory having the same type of values).				*
 *									*
 * writeVector								*
 * -----------								*
 *	Writes a vector to a stream.					*
 *									*
 * checkVector								*
 * -----------								*
 *	Checks if a vector has values of a known type.			*
 *									*
 * zeroVector								*
 * ----------								*
 *	Sets the values in a vector to zero.				*
 *									*
 * allocateVector							*
 * --------------							*
 *	Returns a pointer to a vector of the specified type		*
 *	and number of values, with its values set to zero.		*
 *									*
 * deallocateVector							*
 * ----------------							*
 *	Deallocates memory for a vector.				*
 *									*
 * getVectorValue							*
 * --------------							*
 *	Retrieves a value from a vector, and				*
 *	stores the value in a scalar.					*
 *									*
 * setVectorValue							*
 * --------------							*
 *	Converts a scalar to the type of a vector, and			*
 *	stores the resulting value in the vector.			*
 *									*
 * getVectorValueType							*
 * ------------------							*
 *	Returns the type of values in a vector.				*
 *									*
 * setVectorValueType							*
 * ------------------							*
 *	Sets the type of values in a vector.				*
 *									*
 * getVectorValueCount							*
 * -------------------							*
 *	Returns the number of values in a vector.			*
 *									*
 * setVectorValueCount							*
 * -------------------							*
 *	Sets the number of values in a vector.				*
 *									*
 * getVectorValuePtr							*
 * -----------------							*
 *	Returns a void pointer to a value in a vector.			*
 *									*
 * setVectorValuePtr							*
 * -----------------							*
 *	Sets the pointer to the values in a vector.			*
 *									*
 ************************************************************************/



#include <stdio.h>
#include <memory.h>

/*
 * "old"
 *
#include <malloc.h>
#include <values.h>
 *
 */

/*
 * for MacOS X
 */
#include <stdlib.h>
#include <limits.h>

#include "cfiError.h"
#include "cfiTypes.h"
#include "cfiVector.h"



static char l_msg[MSG_LENGTH];		/* for compound messages */



/************************************************************************
 *									*
 * vectorDotVector							*
 * ---------------							*
 *	Calculates the inner product of two vectors			*
 *	(accelerated for floating-point operands			*
 *	having the same type of values).				*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
vectorDotVector
(
	const cfiVector *const vector1,		/* one vector              */
	const cfiVector *const vector2,		/* the other vector        */
	cfiScalar       *const product_ptr	/* their inner product     */
)

{
	cfiErrorStatus errorStatus;		/* the return value        */

	cfiCounter  valueCount;			/* number of vector values */
	cfiDataType value1Type;			/* types of vector values  */
	cfiDataType value2Type;

	void         *void1_ptr;		/* for accelerated loop    */
	void         *void2_ptr;
	const float  *float1_ptr;
	const float  *float2_ptr;
	const float  *end_float1;
	const double *double1_ptr;
	const double *double2_ptr;
	const double *end_double1;

	cfiScalar  value1;			/* for unaccelerated loop  */
	cfiScalar  value2;
	cfiCounter i;

	DEBUG_MSG(cfiTraceDebug, "calculating inner product of values");

	void1_ptr = getVectorValuePtr(vector1, 0);
	void2_ptr = getVectorValuePtr(vector2, 0);
	value1Type = getVectorValueType(vector1);
	value2Type = getVectorValueType(vector2);
	valueCount = getVectorValueCount(vector1);

	if ( (g_debugStatus & cfiArgDebug) )
	{
		if ( ( void1_ptr == (void *)NULL )
		||   ( void2_ptr == (void *)NULL )
		||   ( value1Type == unknownType )
		||   ( value2Type == unknownType )
		||   ( valueCount != getVectorValueCount(vector2) )
		||   ( valueCount <= 0 )
		||   ( product_ptr == (cfiScalar *)NULL ) )
		{
			FAILURE_MSG("cannot calculate inner product of values");
			return (cfiFailure);
		}
	}

	/*
	 * accelerate calculation for floating-point vectors
	 */

	(*product_ptr) = 0.0;

	if ( ( value1Type == floatType )
	&&   ( value2Type == floatType ) )
	{
		DEBUG_MSG(cfiTraceDebug, "(accelerating float-float inner product)");

		float1_ptr = (float *)void1_ptr;
		float2_ptr = (float *)void2_ptr;

		end_float1 = float1_ptr + valueCount;
		while ( float1_ptr < end_float1 )
		{
			(*product_ptr) += (*float1_ptr++) * (*float2_ptr++);
		}

		return (cfiSuccess);
	}

	if ( ( value1Type == doubleType )
	&&   ( value2Type == doubleType ) )
	{
		DEBUG_MSG(cfiTraceDebug, "(accelerating double-double inner product)");
		double1_ptr = (double *)void1_ptr;
		double2_ptr = (double *)void2_ptr;

		end_double1 = double1_ptr + valueCount;
		while ( double1_ptr < end_double1 )
		{
			(*product_ptr) += (*double1_ptr++) * (*double2_ptr++);
		}

		return (cfiSuccess);
	}

	/*
	 * do an unaccelerated calculation
	 */

	errorStatus = cfiSuccess;

	for ( i = 0; i < valueCount; i++ )
	{
		errorStatus |= getVectorValue(vector1, i, &value1);
		errorStatus |= getVectorValue(vector2, i, &value2);

		(*product_ptr) += value1 * value2;
	}

	if ( (errorStatus & cfiFailure) )
	{
		FAILURE_MSG("cannot calculate inner product of values");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of vectorDotVector() */



/************************************************************************
 *									*
 * vectorOpVector							*
 * --------------							*
 *	Performs pointwise arithmetic operations on two vectors		*
 *	(accelerated for floating-point operands and result		*
 *	having the same type of values).				*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
vectorOpVector
(
	const cfiVector *const vector1,		/* one vector              */
	const cfiOp      op,			/* the operation           */
	const cfiVector *const vector2,		/* the other vector        */
	cfiVector       *const resultVector	/* holds the result        */
)

{
	cfiErrorStatus errorStatus;		/* the return value        */

	cfiCounter  valueCount;			/* number of vector values */
	cfiDataType value1Type;			/* types of vector values  */
	cfiDataType value2Type;
	cfiDataType resultType;

	void         *void1_ptr;		/* for accelerated loop    */
	void         *void2_ptr;
	void         *voidResult_ptr;
	const float  *float1_ptr;
	const float  *float2_ptr;
	const float  *end_float1;
	float        *floatResult_ptr;
	const double *double1_ptr;
	const double *double2_ptr;
	const double *end_double1;
	double       *doubleResult_ptr;

	cfiScalar  value1;			/* for unaccelerated loop  */
	cfiScalar  value2;
	cfiCounter i;

	void1_ptr = getVectorValuePtr(vector1, 0);
	void2_ptr = getVectorValuePtr(vector2, 0);
	voidResult_ptr = getVectorValuePtr(resultVector, 0);
	value1Type = getVectorValueType(vector1);
	value2Type = getVectorValueType(vector2);
	resultType = getVectorValueType(resultVector);
	valueCount = getVectorValueCount(vector1);

	if ( (g_debugStatus & cfiArgDebug) )
	{
		if ( ( void1_ptr == (void *)NULL )
		||   ( void2_ptr == (void *)NULL )
		||   ( voidResult_ptr == (void *)NULL )
		||   ( value1Type == unknownType )
		||   ( value2Type == unknownType )
		||   ( resultType == unknownType )
		||   ( valueCount != getVectorValueCount(vector2) )
		||   ( valueCount != getVectorValueCount(resultVector) )
		||   ( valueCount <= 0 ) )
		{
			FAILURE_MSG("cannot perform operation on values");
			return (cfiFailure);
		}
	}

	/*
 	 * accelerate calculation for floating-point operands
	 * and result having the same type of values
	 */

	if ( ( value1Type == floatType )
	&&   ( value2Type == floatType )
	&&   ( resultType == floatType ) )
	{
		float1_ptr = (float *)void1_ptr;
		float2_ptr = (float *)void2_ptr;
		end_float1 = float1_ptr + valueCount;
		floatResult_ptr = (float *)voidResult_ptr;

		switch (op)
		{
		  case (addOp):

			DEBUG_MSG(cfiTraceDebug, "adding values");
			DEBUG_MSG(cfiTraceDebug, "(accelerating float-float addition)");

			while ( float1_ptr < end_float1 )
			{
				(*floatResult_ptr++) = (*float1_ptr++) + (*float2_ptr++);
			}

			return (cfiSuccess);

		  case (subtractOp):

			DEBUG_MSG(cfiTraceDebug, "subtracting values");
			DEBUG_MSG(cfiTraceDebug, "(accelerating float-float subtraction)");

			while ( float1_ptr < end_float1 )
			{
				(*floatResult_ptr++) = (*float1_ptr++) - (*float2_ptr++);
			}

			return (cfiSuccess);

		  case (multiplyOp):

			DEBUG_MSG(cfiTraceDebug, "multiplying values");
			DEBUG_MSG(cfiTraceDebug, "(accelerating float-float multiplication)");

			while ( float1_ptr < end_float1 )
			{
				(*floatResult_ptr++) = (*float1_ptr++) * (*float2_ptr++);
			}

			return (cfiSuccess);

		  case (divideOp):

			DEBUG_MSG(cfiTraceDebug, "dividing values");
			DEBUG_MSG(cfiTraceDebug, "(accelerating float-float division)");

			while ( float1_ptr < end_float1 )
			{
				if ( *float2_ptr != 0.0 )
				{
					(*floatResult_ptr++) = (*float1_ptr++) / (*float2_ptr++);
				}
				else
				{
					FAILURE_MSG("cannot perform operation on values");
					return (cfiFailure);
				}
			}

			return (cfiSuccess);

		  default:

			FAILURE_MSG("cannot perform operation on values");
			return (cfiFailure);

		} /* end of switch (op) */

	} /* end of if ( ( value1Type == floatType ) ... ) statement */

	if ( ( value1Type == doubleType )
	&&   ( value2Type == doubleType )
	&&   ( resultType == doubleType ) )
	{
		double1_ptr = (double *)void1_ptr;
		double2_ptr = (double *)void2_ptr;
		end_double1 = double1_ptr + valueCount;
		doubleResult_ptr = (double *)voidResult_ptr;

		switch (op)
		{
		  case (addOp):

			DEBUG_MSG(cfiTraceDebug, "adding values");
			DEBUG_MSG(cfiTraceDebug, "(accelerating double-double addition)");

			while ( double1_ptr < end_double1 )
			{
				(*doubleResult_ptr++) = (*double1_ptr++) + (*double2_ptr++);
			}

			return (cfiSuccess);

		  case (subtractOp):

			DEBUG_MSG(cfiTraceDebug, "subtracting values");
			DEBUG_MSG(cfiTraceDebug, "(accelerating double-double subtraction)");

			while ( double1_ptr < end_double1 )
			{
				(*doubleResult_ptr++) = (*double1_ptr++) - (*double2_ptr++);
			}

			return (cfiSuccess);

		  case (multiplyOp):

			DEBUG_MSG(cfiTraceDebug, "multiplying values");
			DEBUG_MSG(cfiTraceDebug, "(accelerating double-double multiplication)");

			while ( double1_ptr < end_double1 )
			{
				(*doubleResult_ptr++) = (*double1_ptr++) * (*double2_ptr++);
			}

			return (cfiSuccess);

		  case (divideOp):

			DEBUG_MSG(cfiTraceDebug, "dividing values");
			DEBUG_MSG(cfiTraceDebug, "(accelerating double-double division)");

			while ( double1_ptr < end_double1 )
			{
				if ( *double2_ptr != 0.0 )
				{
					(*doubleResult_ptr++) = (*double1_ptr++) / (*double2_ptr++);
				}
				else
				{
					FAILURE_MSG("cannot perform operation on values");
					return (cfiFailure);
				}
			}

			return (cfiSuccess);

		  default:

			FAILURE_MSG("cannot perform operation on values");
			return (cfiFailure);

		} /* end of switch (op) */

	} /* end of if ( ( value1Type == doubleType ) ... ) statement */

	/*
	 * do an unaccelerated calculation
	 */

	errorStatus = cfiSuccess;

	switch (op)
	{
	  case (addOp):

		DEBUG_MSG(cfiTraceDebug, "adding values");

		for ( i = 0; i < valueCount; i++ )
		{
			errorStatus |= getVectorValue(vector1, i, &value1);
			errorStatus |= getVectorValue(vector2, i, &value2);
			errorStatus |= setVectorValue(value1 + value2, resultVector, i);
		}

		break;

	  case (subtractOp):

		DEBUG_MSG(cfiTraceDebug, "subtracting values");

		for ( i = 0; i < valueCount; i++ )
		{
			errorStatus |= getVectorValue(vector1, i, &value1);
			errorStatus |= getVectorValue(vector2, i, &value2);
			errorStatus |= setVectorValue(value1 - value2, resultVector, i);
		}

		break;

	  case (multiplyOp):

		DEBUG_MSG(cfiTraceDebug, "multiplying values");

		for ( i = 0; i < valueCount; i++ )
		{
			errorStatus |= getVectorValue(vector1, i, &value1);
			errorStatus |= getVectorValue(vector2, i, &value2);
			errorStatus |= setVectorValue(value1 * value2, resultVector, i);
		}

		break;

	  case (divideOp):

		DEBUG_MSG(cfiTraceDebug, "dividing values");

		for ( i = 0; i < valueCount; i++ )
		{
			errorStatus |= getVectorValue(vector1, i, &value1);
			errorStatus |= getVectorValue(vector2, i, &value2);

			if ( value2 != 0.0 )
			{
				errorStatus |= setVectorValue(value1 / value2, resultVector, i);
			}
			else
			{
				FAILURE_MSG("cannot perform operation on values");
				return (cfiFailure);
			}
		}

		break;

	  default:

		FAILURE_MSG("cannot perform operation on values");
		return (cfiFailure);

	} /* end of switch (op) */

	if ( (errorStatus & cfiFailure) )
	{
		FAILURE_MSG("cannot perform operation on values");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of vectorOpVector() */



/************************************************************************
 *									*
 * scaleVector								*
 * -----------								*
 *	Multiplies a vector by a scalar					*
 *	(accelerated for a floating-point vector).			*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
scaleVector
(
	const cfiScalar  scalar,	/* the scalar              */
	cfiVector       *const vector	/* the vector of interest  */
)

{
	cfiErrorStatus errorStatus;	/* the return value        */

	cfiCounter  valueCount;		/* number of vector values */
	cfiDataType valueType;		/* type of vector values   */

	void   *void_ptr;		/* for accelerated loop    */
	float  *float_ptr;
	float  *end_float;
	double *double_ptr;
	double *end_double;

	cfiScalar  value;		/* for unaccelerated loop  */
	cfiCounter i;

	DEBUG_MSG(cfiTraceDebug, "scaling values");

	void_ptr = getVectorValuePtr(vector, 0);
	valueType = getVectorValueType(vector);
	valueCount = getVectorValueCount(vector);

	if ( (g_debugStatus & cfiArgDebug) )
	{
		if ( ( void_ptr == (void *)NULL )
		||   ( valueType == unknownType )
		||   ( valueCount <= 0 ) )
		{
			FAILURE_MSG("cannot scale values");
			return (cfiFailure);
		}
	}

	/*
	 * accelerate calculation for floating-point vector
	 */

	switch (valueType)
	{
	  case floatType:

		DEBUG_MSG(cfiTraceDebug, "(accelerating float scaling)");

		float_ptr = (float *)void_ptr;
		end_float = float_ptr + valueCount;
		while ( float_ptr < end_float )
		{
			(*float_ptr++) *= scalar;
		}

		return (cfiSuccess);

	  case doubleType:

		DEBUG_MSG(cfiTraceDebug, "(accelerating double scaling)");

		double_ptr = (double *)void_ptr;
		end_double = double_ptr + valueCount;
		while ( double_ptr < end_double )
		{
			(*double_ptr++) *= scalar;
		}

		return (cfiSuccess);

	  default:

		break;		/* do an unaccelerated calculation */

	} /* end of switch (valueType) */

	/*
	 * do an unaccelerated calculation
	 */

	errorStatus = cfiSuccess;

	for ( i = 0; i < valueCount; i++ )
	{
		errorStatus |= getVectorValue(vector, i, &value);
		errorStatus |= setVectorValue(scalar * value, vector, i);
	}

	if ( (errorStatus & cfiFailure) )
	{
		FAILURE_MSG("cannot scale values");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of scaleVector() */



/************************************************************************
 *									*
 * copyVector								*
 * ----------								*
 *	Copies values from a vector of arbitrary			*
 *	type to another vector of arbitrary type			*
 *	(accelerated for vectors in non-overlapping			*
 *	memory having the same type of values).				*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
copyVector
(
	const cfiVector *const fromVector,	/* the source vector       */
	cfiVector       *const toVector		/* the destination vector  */
)

{
	cfiErrorStatus errorStatus;		/* the return value        */

	cfiCounter  valueCount;			/* number of vector values */
	cfiDataType fromType;			/* types of vector values  */
	cfiDataType toType;

	void *from_ptr;				/* for accelerated copy    */
	void *end_from;
	void *to_ptr;
	void *end_to;

	cfiScalar  value;			/* for unaccelerated copy  */
	cfiCounter i;

	DEBUG_MSG(cfiTraceDebug, "copying values");

	valueCount = getVectorValueCount(fromVector);
	from_ptr = getVectorValuePtr(fromVector, 0);
	end_from = getVectorValuePtr(fromVector, valueCount);
	to_ptr = getVectorValuePtr(toVector, 0);
	end_to = getVectorValuePtr(toVector, valueCount);
	fromType = getVectorValueType(fromVector);
	toType = getVectorValueType(toVector);

	if ( (g_debugStatus & cfiArgDebug) )
	{
		if ( ( valueCount != getVectorValueCount(toVector) )
		||   ( valueCount <= 0 )
		||   ( from_ptr == (void *)NULL )
		||   ( end_from == (void *)NULL )
		||   ( to_ptr == (void *)NULL )
		||   ( end_to == (void *)NULL )
		||   ( fromType == unknownType )
		||   ( toType == unknownType ) )
		{
			FAILURE_MSG("cannot copy values");
			return (cfiFailure);
		}
	}

	errorStatus = cfiSuccess;

	if ( ( ( to_ptr >= from_ptr ) && ( to_ptr < end_from ) )
	||   ( ( from_ptr >= to_ptr ) && ( from_ptr < end_to ) ) )
	{
		WARNING_MSG("copying values in overlapping memory");
		errorStatus |= cfiWarning;
	}
	else if ( toType == fromType )
	{
		/*
		 * accelerate copying of non-overlapping
		 * values of the same type
		 */

		DEBUG_MSG(cfiTraceDebug, "(accelerating copying of values of same type)");
		(void)memcpy(to_ptr, from_ptr, ((char *)end_from) - ((char *)from_ptr));
		return (cfiSuccess);
	}

	/*
	 * do an unaccelerated copy
	 */

	for ( i = 0; i < valueCount; i++ )
	{
		errorStatus |= getVectorValue(fromVector, i, &value);
		errorStatus |= setVectorValue(value, toVector, i);
	}

	if ( (errorStatus & cfiFailure) )
	{
		FAILURE_MSG("cannot copy values");
		return (cfiFailure);
	}

	return (errorStatus);

} /* end of copyVector() */



/************************************************************************
 *									*
 * writeVector								*
 * -----------								*
 *	Writes a vector to a stream.					*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
writeVector
(
	const cfiVector  *const vector,		/* the vector of interest    */
	const cfiCounter  valuesPerLine, 	/* number of values per line */
	FILE             *const stream		/* the output stream         */
)

{
	cfiCounter valueCount;			/* number of vector values   */
	cfiCounter i;				/* loop control              */

	void                *void_ptr;		/* loop control              */
	const unsigned char *byte_ptr;
	const int         *short_ptr;
	const long          *long_ptr;
	const float         *float_ptr;
	const double        *double_ptr;

	(void)sprintf(l_msg, "writing %d value(s) per line", valuesPerLine);
	DEBUG_MSG(cfiInOutDebug, l_msg);

	if ( valuesPerLine < 0 )
	{
		FAILURE_MSG("cannot write values");
		return (cfiFailure);
	}

	if ( vector == (cfiVector *)NULL )
	{
		if ( fprintf(stream, "null vector\n") < 0 )
		{
			FAILURE_MSG("cannot write values");
			return (cfiFailure);
		}

		return (cfiSuccess);
	}

	void_ptr = getVectorValuePtr(vector, 0);
	valueCount = getVectorValueCount(vector);

	switch (getVectorValueType(vector))
	{
	  case byteType:

		(void)fprintf(stream, "%d byteType [", valueCount);

		if ( ( valuesPerLine > 0 )
		&&   ( void_ptr != (void *)NULL ) )
		{
			byte_ptr = (unsigned char *)void_ptr;
			for ( i = 0; i < valueCount; i++ )
			{
				if ( (i % valuesPerLine) == 0 )
				{
					(void)fprintf(stream, "\n");
				}
				(void)fprintf(stream, " %u", *byte_ptr++);
			}
		}

		(void)fprintf(stream, "\n]\n");
		break;

	  case shortType:

		(void)fprintf(stream, "%d shortType [", valueCount);

		if ( ( valuesPerLine > 0 )
		&&   ( void_ptr != (void *)NULL ) )
		{
			short_ptr = (int *)void_ptr;
			for ( i = 0; i < valueCount; i++ )
			{
				if ( (i % valuesPerLine) == 0 )
				{
					(void)fprintf(stream, "\n");
				}
				(void)fprintf(stream, " %hd", *short_ptr++);
			}
		}

		(void)fprintf(stream, "\n]\n");
		break;

	  case longType:

		(void)fprintf(stream, "%d longType [", valueCount);

		if ( ( valuesPerLine > 0 )
		&&   ( void_ptr != (void *)NULL ) )
		{
			long_ptr = (long *)void_ptr;
			for ( i = 0; i < valueCount; i++ )
			{
				if ( (i % valuesPerLine) == 0 )
				{
					(void)fprintf(stream, "\n");
				}
				(void)fprintf(stream, " %ld", *long_ptr++);
			}
		}

		(void)fprintf(stream, "\n]\n");
		break;

	  case floatType:

		(void)fprintf(stream, "%d floatType [", valueCount);

		if ( ( valuesPerLine > 0 )
		&&   ( void_ptr != (void *)NULL ) )
		{
			float_ptr = (float *)void_ptr;
			for ( i = 0; i < valueCount; i++ )
			{
				if ( (i % valuesPerLine) == 0 )
				{
					(void)fprintf(stream, "\n");
				}
				(void)fprintf(stream, " %G", *float_ptr++);
			}
		}

		(void)fprintf(stream, "\n]\n");
		break;

	  case doubleType:

		(void)fprintf(stream, "%d doubleType [", valueCount);

		if ( ( valuesPerLine > 0 )
		&&   ( void_ptr != (void *)NULL ) )
		{
			double_ptr = (double *)void_ptr;
			for ( i = 0; i < valueCount; i++ )
			{
				if ( (i % valuesPerLine) == 0 )
				{
					(void)fprintf(stream, "\n");
				}
				(void)fprintf(stream, " %G", *double_ptr++);
			}
		}

		(void)fprintf(stream, "\n]\n");
		break;

	  default:

		(void)fprintf(stream, "%d unknownType [", valueCount);
		(void)fprintf(stream, "\n]\n");
		break;

	} /* end of switch (getVectorValueType(vector)) */

	if ( ferror(stream) )
	{
		FAILURE_MSG("cannot write values");
		return (cfiFailure);
	}

	return (cfiSuccess);

} /* end of writeVector() */



/************************************************************************
 *									*
 * checkVector								*
 * -----------								*
 *	Checks if a vector has values of a known type.			*
 *									*
 ********************************************************************bwr*/

cfiVectorStatus
checkVector
(
	const cfiVector *const vector	/* the vector of interest */
)

{
	cfiVectorStatus vectorStatus;	/* the return value       */

	DEBUG_MSG(cfiTraceDebug, "checking values");

	vectorStatus = emptyVector;

	if ( ( getVectorValueType(vector) != unknownType )
	&&   ( getVectorValueCount(vector) > 0 )
	&&   ( getVectorValuePtr(vector, 0) != (void *)NULL ) )
	{
		vectorStatus |= vectorHasValues;
	}

	return (vectorStatus);

} /* end of checkVector() */



/************************************************************************
 *									*
 * zeroVector								*
 * ----------								*
 *	Sets the values in a vector to zero.				*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
zeroVector
(
	cfiVector *const vector		/* the vector of interest  */
)

{
	cfiCounter valueCount;		/* number of vector values */

	void          *void_ptr;	/* loop control            */
	unsigned char *byte_ptr;
	unsigned char *end_byte;
	int         *short_ptr;
	int         *end_short;
	long          *long_ptr;
	long          *end_long;
	float         *float_ptr;
	float         *end_float;
	double        *double_ptr;
	double        *end_double;

	DEBUG_MSG(cfiTraceDebug, "setting values to zero");

	void_ptr = getVectorValuePtr(vector, 0);
	valueCount = getVectorValueCount(vector);

	if ( (g_debugStatus & cfiArgDebug) )
	{
		if ( ( void_ptr == (void *)NULL )
		||   ( valueCount <= 0 ) )
		{
			FAILURE_MSG("cannot set values to zero");
			return (cfiFailure);
		}
	}

	switch (getVectorValueType(vector))
	{
	  case byteType:

		byte_ptr = (unsigned char *)void_ptr;
		end_byte = byte_ptr + valueCount;
		while ( byte_ptr < end_byte )
		{
			(*byte_ptr++) = 0;
		}

		return (cfiSuccess);

	  case shortType:

		short_ptr = (int *)void_ptr;
		end_short = short_ptr + valueCount;
		while ( short_ptr < end_short )
		{
			(*short_ptr++) = 0;
		}

		return (cfiSuccess);

	  case longType:

		long_ptr = (long *)void_ptr;
		end_long = long_ptr + valueCount;
		while ( long_ptr < end_long )
		{
			(*long_ptr++) = 0;
		}

		return (cfiSuccess);

	  case floatType:

		float_ptr = (float *)void_ptr;
		end_float = float_ptr + valueCount;
		while ( float_ptr < end_float )
		{
			(*float_ptr++) = 0.0;
		}

		return (cfiSuccess);

	  case doubleType:

		double_ptr = (double *)void_ptr;
		end_double = double_ptr + valueCount;
		while ( double_ptr < end_double )
		{
			(*double_ptr++) = 0.0;
		}

		return (cfiSuccess);

	  default:

		break;

	} /* end of switch (getVectorValueType(vector)) */

	FAILURE_MSG("cannot set values to zero");
	return (cfiFailure);

} /* end of zeroVector() */



/************************************************************************
 *									*
 * allocateVector							*
 * --------------							*
 *	Returns a pointer to a vector of the specified type		*
 *	and number of values, with its values set to zero.		*
 *									*
 ********************************************************************bwr*/

cfiVector *
allocateVector
(
	const cfiDataType valueType,	/* type of vector values        */
	const cfiCounter  valueCount	/* number of vector values      */
)

{
	cfiVector *vector;		/* the return value             */

	cfiCounter memorySize;		/* amount of memory to allocate */

	(void)sprintf(l_msg, "allocating %d value(s)", valueCount);
	DEBUG_MSG(cfiMemoryDebug, l_msg);

	if ( valueCount < 0 )
	{
		FAILURE_MSG("cannot allocate values");
		return ((cfiVector *)NULL);
	}

	switch (valueType)
	{
	  case byteType:

		DEBUG_MSG(cfiMemoryDebug, "...of type byte");
		memorySize = valueCount * sizeof(unsigned char);
		break;

	  case shortType:

		DEBUG_MSG(cfiMemoryDebug, "...of type int");
		memorySize = valueCount * sizeof(int);
		break;

	  case longType:

		DEBUG_MSG(cfiMemoryDebug, "...of type long");
		memorySize = valueCount * sizeof(long);
		break;

	  case floatType:

		DEBUG_MSG(cfiMemoryDebug, "...of type float");
		memorySize = valueCount * sizeof(float);
		break;

	  case doubleType:

		DEBUG_MSG(cfiMemoryDebug, "...of type double");
		memorySize = valueCount * sizeof(double);
		break;

	  case unknownType:

		DEBUG_MSG(cfiMemoryDebug, "...of type unknown");
		memorySize = 0;
		break;

	  default:

		DEBUG_MSG(cfiMemoryDebug, "...of ill-defined type");
		FAILURE_MSG("cannot allocate values");
		return ((cfiVector *)NULL);

	} /* end of switch (valueType) */

/* LINTED */
	if ( (vector = (cfiVector *)malloc(sizeof(cfiVector)))
		== (cfiVector *)NULL )
	{
		FAILURE_MSG("cannot allocate values");
		return ((cfiVector *)NULL);
	}

	vector->type = valueType;
	vector->count = valueCount;
	vector->value = (void *)NULL;

	if ( memorySize > 0 )
	{
		if ( ( ( vector->value = (void *)malloc( ( size_t ) memorySize ) )
			== (void *) NULL )
		||   ( zeroVector(vector) != cfiSuccess ) )
		{
			FAILURE_MSG("cannot allocate values");
			deallocateVector(&vector);
			return ((cfiVector *)NULL);
		}
	}

	return (vector);

} /* end of allocateVector() */



/************************************************************************
 *									*
 * deallocateVector							*
 * ----------------							*
 *	Deallocates memory for a vector.				*
 *									*
 ********************************************************************bwr*/

void
deallocateVector
(
	cfiVector **const vector_ptr	/* the vector of interest */
)

{
	DEBUG_MSG(cfiMemoryDebug, "deallocating values");

	if ( (*vector_ptr) == (cfiVector *)NULL )
	{
		return;
	}

	if ( (*vector_ptr)->value != (void *)NULL )
	{
		(void)free((*vector_ptr)->value);
	}

	(void)free(*vector_ptr);
	(*vector_ptr) = (cfiVector *)NULL;

} /* end of deallocateVector() */



/************************************************************************
 *									*
 * getVectorValue							*
 * --------------							*
 *	Retrieves a value from a vector, and				*
 *	stores the value in a scalar.					*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
getVectorValue
(
	const cfiVector  *const vector,		/* the source vector          */
	const cfiCounter  index,		/* index of value of interest */
	cfiScalar        *const value_ptr	/* destination of the value   */
)

{
	void *void_ptr;				/* pointer to vector values   */

	if ( (g_debugStatus & cfiArgDebug) )
	{
		if ( ( vector == (cfiVector *)NULL )
		||   ( (void_ptr = vector->value)
			== (void *)NULL )
		||   ( index < 0 )
		||   ( index >= vector->count )
		||   ( value_ptr == (cfiScalar *)NULL ) )
		{
			FAILURE_MSG("cannot get value");
			return (cfiFailure);
		}
	}
	else
	{
		void_ptr = vector->value;
	}

	switch (vector->type)
	{
	  case byteType:

		(*value_ptr) = (cfiScalar)(*(((unsigned char *)void_ptr) + index));
		return (cfiSuccess);

	  case shortType:

		(*value_ptr) = (cfiScalar)(*(((int *)void_ptr) + index));
		return (cfiSuccess);

	  case longType:

		(*value_ptr) = (cfiScalar)(*(((long *)void_ptr) + index));
		return (cfiSuccess);

	  case floatType:

		(*value_ptr) = (cfiScalar)(*(((float *)void_ptr) + index));
		return (cfiSuccess);

	  case doubleType:

		(*value_ptr) = (cfiScalar)(*(((double *)void_ptr) + index));
		return (cfiSuccess);

	  default:

		break;

	} /* end of switch (vector->type) */

	FAILURE_MSG("cannot get value");
	return (cfiFailure);

} /* end of getVectorValue() */



/************************************************************************
 *									*
 * setVectorValue							*
 * --------------							*
 *	Converts a scalar to the type of a vector, and			*
 *	stores the resulting value in the vector.			*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
setVectorValue
(
	cfiScalar         value,	/* value to be stored in vector */
	cfiVector        *const vector,	/* the destination vector       */
	const cfiCounter  index		/* index of value of interest   */
)

{
	void *void_ptr;			/* pointer to vector values     */

	if ( (g_debugStatus & cfiArgDebug) )
	{
		if ( ( vector == (cfiVector *)NULL )
		||   ( (void_ptr = vector->value)
			== (void *)NULL )
		||   ( index < 0 )
		||   ( index >= vector->count ) )
		{
			FAILURE_MSG("cannot set value");
			return (cfiFailure);
		}
	}
	else
	{
		void_ptr = vector->value;
	}

	switch (vector->type)
	{
	  case byteType:

		if ( value > 255.0 )
		{
			DEBUG_MSG(cfiMaxDebug, "clipping value of type byte");
			value = 255.0;
		}
		else if ( value < 0.0 )
		{
			DEBUG_MSG(cfiMaxDebug, "clipping value of type byte");
			value = 0.0;
		}
		else
		{
			value += 0.5;
		}

		*(((unsigned char *)void_ptr) + index) = (unsigned char)value;
		return (cfiSuccess);

	  case shortType:

		if ( value > (cfiScalar)MAXSHORT )
		{
			DEBUG_MSG(cfiMaxDebug, "clipping value of type int");
			value = (cfiScalar)MAXSHORT;
		}
		else if ( value < (cfiScalar)(-MAXSHORT) )
		{
			DEBUG_MSG(cfiMaxDebug, "clipping value of type int");
			value = (cfiScalar)(-MAXSHORT);
		}
		else if ( value >= 0.0 )
		{
			value += 0.5;
		}
		else
		{
			value -= 0.5;
		}

		*(((int *)void_ptr) + index) = (int)value;
		return (cfiSuccess);

	  case longType:

/* LINTED */
		if ( value > (cfiScalar)MAXLONG )
		{
			DEBUG_MSG(cfiMaxDebug, "clipping value of type long");
/* LINTED */
			value = (cfiScalar)MAXLONG;
		}
/* LINTED */
		else if ( value < (cfiScalar)(-MAXLONG) )
		{
			DEBUG_MSG(cfiMaxDebug, "clipping value of type long");
/* LINTED */
			value = (cfiScalar)(-MAXLONG);
		}
		else if ( value >= 0.0 )
		{
			value += 0.5;
		}
		else
		{
			value -= 0.5;
		}

		*(((long *)void_ptr) + index) = (long)value;
		return (cfiSuccess);

	  case floatType:

		if ( value > (cfiScalar)MAXFLOAT )
		{
			DEBUG_MSG(cfiMaxDebug, "clipping value of type float");
			value = (cfiScalar)MAXFLOAT;
		}
		else if ( value < (cfiScalar)(-MAXFLOAT) )
		{
			DEBUG_MSG(cfiMaxDebug, "clipping value of type float");
			value = (cfiScalar)(-MAXFLOAT);
		}

		*(((float *)void_ptr) + index) = (float)value;
		return (cfiSuccess);

	  case doubleType:

		*(((double *)void_ptr) + index) = (double)value;
		return (cfiSuccess);

	  default:

		break;

	} /* end of switch (vector->type) */

	FAILURE_MSG("cannot set value");
	return (cfiFailure);

} /* end of setVectorValue() */



/************************************************************************
 *									*
 * getVectorValueType							*
 * ------------------							*
 *	Returns the type of values in a vector.				*
 *									*
 ********************************************************************bwr*/

cfiDataType
getVectorValueType
(
	const cfiVector *const vector	/* the vector of interest */
)

{
	cfiDataType valueType;		/* the return value       */

	if ( vector == (cfiVector *)NULL )
	{
		DEBUG_MSG(cfiReturnDebug, "type of values is unknown");
		return (unknownType);
	}

	valueType = vector->type;

	if ( ( valueType == byteType )
	||   ( valueType == shortType )
	||   ( valueType == longType )
	||   ( valueType == floatType )
	||   ( valueType == doubleType ) )
	{
		return (valueType);
	}

	DEBUG_MSG(cfiReturnDebug, "type of values is unknown");
	return (unknownType);

} /* end of getVectorValueType() */



/************************************************************************
 *									*
 * setVectorValueType							*
 * ------------------							*
 *	Sets the type of values in a vector.				*
 *									*
 * NOTE-To maintain a minimum level of data integrity,			*
 *	the vector must have a NULL value pointer.			*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
setVectorValueType
(
	const cfiDataType  valueType,	/* the desired type of values */
	cfiVector         *const vector	/* the vector of interest     */
)

{
	if ( ( vector == (cfiVector *)NULL )
	||   ( vector->value != (void *)NULL ) )
	{
		FAILURE_MSG("cannot set type of values");
		return (cfiFailure);
	}

	if ( ( valueType == (cfiDataType)byteType )
	||   ( valueType == (cfiDataType)shortType )
	||   ( valueType == (cfiDataType)longType )
	||   ( valueType == (cfiDataType)floatType )
	||   ( valueType == (cfiDataType)doubleType )
	||   ( valueType == (cfiDataType)unknownType ) )
	{
		vector->type = valueType;
		return (cfiSuccess);
	}

	FAILURE_MSG("cannot set type of values");
	return (cfiFailure);

} /* end of setVectorValueType() */



/************************************************************************
 *									*
 * getVectorValueCount							*
 * -------------------							*
 *	Returns the number of values in a vector.			*
 *									*
 ********************************************************************bwr*/

cfiCounter
getVectorValueCount
(
	const cfiVector *const vector	/* the vector of interest */
)

{
	if ( ( vector == (cfiVector *)NULL )
	||   ( vector->count <= 0 ) )
	{
		DEBUG_MSG(cfiReturnDebug, "number of values is ill-defined");
		return (0);
	}

	return (vector->count);

} /* end of getVectorValueCount() */



/************************************************************************
 *									*
 * setVectorValueCount							*
 * -------------------							*
 *	Sets the number of values in a vector.				*
 *									*
 * NOTE-To maintain a minimum level of data integrity,			*
 *	the vector must have a NULL value pointer.			*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
setVectorValueCount
(
	const cfiCounter  valueCount,	/* the desired number of values */
	cfiVector        *const vector	/* the vector of interest       */
)

{
	if ( ( vector == (cfiVector *)NULL )
	||   ( vector->value != (void *)NULL )
	||   ( valueCount < 0 ) )
	{

		FAILURE_MSG("cannot set number of values");
		return (cfiFailure);
	}

	vector->count = valueCount;
	return (cfiSuccess);

} /* end of setVectorValueCount() */



/************************************************************************
 *									*
 * getVectorValuePtr							*
 * -----------------							*
 *	Returns a void pointer to a value in a vector.			*
 *									*
 ********************************************************************bwr*/

void *
getVectorValuePtr
(
	const cfiVector  *const vector,	/* the vector of interest     */
	const cfiCounter  index		/* index of value of interest */
)

{
	if ( ( vector == (cfiVector *)NULL )
	||   ( vector->value == (void *)NULL )
/* LINTED */
	||   ( index == UNDEFINED_COUNT ) )
	{
		DEBUG_MSG(cfiReturnDebug, "cannot get pointer to value");
		return ((void *)NULL);
	}

	if ( index == 0 )
	{
		return (vector->value);
	}

	switch (vector->type)
	{
	  case byteType:

		return ((void *)(((unsigned char *)vector->value) + index));

	  case shortType:

		return ((void *)(((int *)vector->value) + index));

	  case longType:

		return ((void *)(((long *)vector->value) + index));

	  case floatType:

		return ((void *)(((float *)vector->value) + index));

	  case doubleType:

		return ((void *)(((double *)vector->value) + index));

	  default:

		break;

	} /* end of switch (vector->type) */

	DEBUG_MSG(cfiReturnDebug, "cannot get pointer to value");
	return ((void *)NULL);

} /* end of getVectorValuePtr() */



/************************************************************************
 *									*
 * setVectorValuePtr							*
 * -----------------							*
 *	Sets the pointer to the values in a vector.			*
 *									*
 ********************************************************************bwr*/

cfiErrorStatus
setVectorValuePtr
(
	void      *const value_ptr,	/* pointer to values of interest */
	cfiVector *const vector		/* the vector of interest        */
)

{
	if ( vector == (cfiVector *)NULL )
	{
		FAILURE_MSG("cannot set pointer to values");
		return (cfiFailure);
	}

	vector->value = value_ptr;
	return (cfiSuccess);

} /* end of setVectorValuePtr() */
