/*
 * SCCS:  @(#)cfiVector.h  1.9  6/7/95  11:33:54
 */

#ifndef cfiVector_HEADER
#define cfiVector_HEADER

#ifdef __cplusplus
extern "C" {
#endif



#include <stdio.h>

#include "cfiError.h"
#include "cfiTypes.h"



typedef unsigned int cfiVectorStatus;

typedef enum cfiVectorFlag
{
	emptyVector     = 0x00,
	vectorHasValues = 0x01

} cfiVectorFlag;



typedef struct cfiVector
{
	cfiDataType  type;	/* type of vector values   */
	cfiCounter   count;	/* number of vector values */
	void        *value;	/* the vector values       */

} cfiVector;



/************************************************************************
 *									*
 * prototypes for functions in cfiVector.c				*
 *									*
 ************************************************************************/

extern cfiErrorStatus	vectorDotVector(
				const cfiVector *const,
				const cfiVector *const,
				cfiScalar *const
			);
extern cfiErrorStatus	vectorOpVector(
				const cfiVector *const,
				const cfiOp,
				const cfiVector *const,
				cfiVector *const
			);
extern cfiErrorStatus	scaleVector(
				const cfiScalar,
				cfiVector *const
			);
extern cfiErrorStatus	copyVector(
				const cfiVector *const,
				cfiVector *const
			);
extern cfiErrorStatus	writeVector(
				const cfiVector *const,
				const cfiCounter,
				FILE *const
			);
extern cfiVectorStatus	checkVector(
				const cfiVector *const
			);
extern cfiErrorStatus	zeroVector(
				cfiVector *const
			);
extern cfiVector *	allocateVector(
				const cfiDataType,
				const cfiCounter
			);
extern void		deallocateVector(
				cfiVector **const
			);
extern cfiErrorStatus	getVectorValue(
				const cfiVector *const,
				const cfiCounter,
				cfiScalar *const
			);
extern cfiErrorStatus	setVectorValue(
				cfiScalar,
				cfiVector *const,
				const cfiCounter
			);
extern cfiDataType	getVectorValueType(
				const cfiVector *const
			);
extern cfiErrorStatus	setVectorValueType(
				const cfiDataType,
				cfiVector *const
			);
extern cfiCounter	getVectorValueCount(
				const cfiVector *const
			);
extern cfiErrorStatus	setVectorValueCount(
				const cfiCounter,
				cfiVector *const
			);
extern void *		getVectorValuePtr(
				const cfiVector *const,
				const cfiCounter
			);
extern cfiErrorStatus	setVectorValuePtr(
				void *const,
				cfiVector *const
			);



#ifdef __cplusplus
}
#endif

#endif /* cfiVector_HEADER */
