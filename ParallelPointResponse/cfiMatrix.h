/*
 * SCCS:  @(#)cfiMatrix.h  1.5  6/7/95  11:33:53
 */

#ifndef cfiMatrix_HEADER
#define cfiMatrix_HEADER

#ifdef __cplusplus
extern "C" {
#endif



#include <stdio.h>

#include "cfiError.h"
#include "cfiTypes.h"
#include "cfiVector.h"



typedef unsigned int cfiMatrixStatus;

typedef enum cfiMatrixFlag
{
	emptyMatrix         = 0x00,
	matrixHasValues     = 0x01,
	matrixHasDimensions = 0x02

} cfiMatrixFlag;



typedef struct cfiMatrix
{
	cfiCounter  rowCount;		/* number of rows in the matrix    */
	cfiCounter  colCount;		/* number of columns in the matrix */
	cfiVector  *matrixValues;	/* stores the matrix values        */

} cfiMatrix;



/************************************************************************
 *									*
 * prototypes for functions in cfiMatrix.c				*
 *									*
 ************************************************************************/

extern cfiErrorStatus	matrixTimesMatrix(
				const cfiMatrix *const,
				const cfiMatrix *const,
				cfiMatrix *const
			);
extern cfiErrorStatus	matrixTimesVector(
				const cfiMatrix *const,
				const cfiVector *const,
				cfiVector *const
			);
extern cfiErrorStatus	vectorTimesMatrix(
				const cfiVector *const,
				const cfiMatrix *const,
				cfiVector *const
			);
extern cfiMatrix *	allocateMatrix(
				const cfiDataType,
				const cfiCounter,
				const cfiCounter
			);
extern void		deallocateMatrix(
				cfiMatrix **const
			);
extern cfiErrorStatus	matrixOpMatrix(
				const cfiMatrix *const,
				const cfiOp,
				const cfiMatrix *const,
				cfiMatrix *const
			);
extern cfiErrorStatus	scaleMatrix(
				const cfiScalar,
				cfiMatrix *const
			);
extern cfiErrorStatus	copyMatrix(
				const cfiMatrix *const,
				cfiMatrix *const
			);
extern cfiErrorStatus	writeMatrix(
				const cfiMatrix *const,
				FILE *const
			);
extern cfiMatrixStatus	checkMatrix(
				const cfiMatrix *const
			);
extern cfiErrorStatus	zeroMatrix(
				cfiMatrix *const
			);
extern cfiErrorStatus	getMatrixValues(
				const cfiMatrix *const,
				const cfiCounter,
				const cfiCounter,
				cfiScalar *
			);
extern cfiErrorStatus	setMatrixValues(
				const cfiScalar *,
				cfiMatrix *const,
				const cfiCounter,
				const cfiCounter
			);
extern cfiDataType	getMatrixValueType(
				const cfiMatrix *const
			);
extern cfiErrorStatus	setMatrixValueType(
				const cfiDataType,
				cfiMatrix *const
			);
extern cfiCounter	getMatrixDimension(
				const cfiMatrix *const,
				const cfiDataAxis
			);
extern cfiErrorStatus	setMatrixDimension(
				const cfiDataAxis,
				const cfiCounter,
				cfiMatrix *const
			);
extern void *		getMatrixValuePtr(
				const cfiMatrix *const,
				const cfiCounter,
				const cfiCounter
			);
extern cfiErrorStatus	setMatrixValuePtr(
				void *const,
				cfiMatrix *const
			);



#ifdef __cplusplus
}
#endif

#endif /* cfiMatrix_HEADER */
