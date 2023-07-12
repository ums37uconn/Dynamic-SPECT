/*
 * SCCS:  @(#)cfiArray.h  1.5  6/7/95  11:33:52
 */

#ifndef cfiArray_HEADER
#define cfiArray_HEADER

#ifdef __cplusplus
extern "C" {
#endif



#include <stdio.h>

#include "cfiError.h"
#include "cfiTypes.h"
#include "cfiVector.h"



typedef unsigned int cfiArrayStatus;

typedef enum cfiArrayFlag
{
	emptyArray         = 0x00,
	arrayHasValues     = 0x01,
	arrayHasDimensions = 0x02

} cfiArrayFlag;



/*
 * define the maximum number of axes for a 2x2x2x...x2 array, such
 * that the total number of values in the array can be counted using
 * the (signed) cfiCounter type defined in cfiTypes.c
 */

/*
#define MAX_AXIS	(BITS(cfiCounter) - 2)
 */
#define MAX_AXIS	(10)



typedef struct cfiArray
{
	cfiCounter  axisCount;			/* number of array axes     */
	cfiCounter  dimensions[MAX_AXIS];	/* #'s of values along axes */
	cfiVector  *arrayValues;		/* stores the array values  */

} cfiArray;



/************************************************************************
 *									*
 * prototypes for functions in cfiArray.c				*
 *									*
 ************************************************************************/

extern cfiErrorStatus	getArrayValues(
				const cfiArray *const,
				const cfiCounter *const,
				cfiScalar *
			);
extern cfiErrorStatus	setArrayValues(
				const cfiScalar *,
				cfiArray *const,
				const cfiCounter *const
			);
extern cfiArray *	allocateArray(
				const cfiDataType,
				const cfiCounter,
				const cfiCounter *const
			);
extern void		deallocateArray(
				cfiArray **const
			);
extern cfiErrorStatus	arrayOpArray(
				const cfiArray *const,
				const cfiOp,
				const cfiArray *const,
				cfiArray *const
			);
extern cfiErrorStatus	scaleArray(
				const cfiScalar,
				cfiArray *const
			);
extern cfiErrorStatus	copyArray(
				const cfiArray *const,
				cfiArray *const
			);
extern cfiErrorStatus	writeArray(
				const cfiArray *const,
				FILE *const
			);
extern cfiArrayStatus	checkArray(
				const cfiArray *const
			);
extern cfiErrorStatus	zeroArray(
				cfiArray *const
			);
extern cfiDataType	getArrayValueType(
				const cfiArray *const
			);
extern cfiErrorStatus	setArrayValueType(
				const cfiDataType,
				cfiArray *const
			);
extern cfiCounter	getArrayAxisCount(
				const cfiArray *const
			);
extern cfiErrorStatus	setArrayAxisCount(
				const cfiCounter,
				cfiArray *const
			);
extern cfiCounter	getArrayDimension(
				const cfiArray *const,
				const cfiDataAxis
			);
extern cfiErrorStatus	setArrayDimension(
				const cfiDataAxis,
				const cfiCounter,
				cfiArray *const
			);
extern void *		getArrayValuePtr(
				const cfiArray *const,
				const cfiCounter *const
			);
extern cfiErrorStatus	setArrayValuePtr(
				void *const,
				cfiArray *const
			);



#ifdef __cplusplus
}
#endif

#endif /* cfiArray_HEADER */
