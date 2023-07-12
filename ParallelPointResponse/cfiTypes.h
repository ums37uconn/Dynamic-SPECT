/*
 * SCCS:  @(#)cfiTypes.h  1.5  4/17/95  10:31:36
 */

#define M_PI 3.14159265358979323846
#define MAXFLOAT ((float)3.40282346638528860e+38)


#ifndef cfiTypes_HEADER
#define cfiTypes_HEADER

#ifdef __cplusplus
extern "C" {
#endif



/*
 * "old"
 *
#include <values.h>
 *
 */

/*
 * for MacOS X
 */
#include <limits.h>
#include <float.h>
#include <math.h>
#define MAXSHORT	SHRT_MAX
#define MAXLONG		LONG_MAX



typedef enum cfiBoolean
{
	cfiFalse,
	cfiTrue

} cfiBoolean;



typedef long   cfiCounter;		/* default fixed-point type     */
typedef double cfiScalar;		/* default floating-point type  */
//typedef float cfiScalar;		/*  */



#define UNDEFINED_COUNT	(~MAXLONG)	/* see cfiCounter typedef above */



typedef enum cfiDataType
{
	unknownType = 0,
	byteType    = 1,
	shortType   = 2,
	longType    = 3,
	floatType   = 4,
	doubleType  = 5,

	counterType = 3,		/* see cfiCounter typedef above */
	scalarType  = 5			/* see cfiScalar typedef above  */

} cfiDataType;



typedef enum cfiOp
{
	addOp,
	subtractOp,
	multiplyOp,
	divideOp

} cfiOp;



typedef enum cfiDataAxis
{
	columnAxis = 0,
	rowAxis    = 1,
	planeAxis  = 2,

	xAxis = 0,
	yAxis = 1,
	zAxis = 2,

	axis0 = 0,
	axis1 = 1,
	axis2 = 2,
	axis3 = 3,
	axis4 = 4,
	axis5 = 5,
	axis6 = 6,
	axis7 = 7,
	axis8 = 8,
	axis9 = 9

} cfiDataAxis;



#ifdef __cplusplus
}
#endif

#endif /* cfiTypes_HEADER */
