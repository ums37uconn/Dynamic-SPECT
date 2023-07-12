/*
 * SCCS:  %Z%%M%  %I%  %G%  %U%
 */

#ifndef cfiMesh_HEADER
#define cfiMesh_HEADER

#ifdef __cplusplus
extern "C" {
#endif



#include <stdio.h>

#include "cfiError.h"
#include "cfiTypes.h"
#include "cfiMatrix.h"
#include "cfiArray.h"



typedef unsigned int cfiMeshStatus;

typedef enum cfiMeshFlag
{
	emptyMesh         = 0x00,
	meshHasValues     = 0x01,
	meshHasDimensions = 0x02,
	meshHasCoords     = 0x04

} cfiMeshFlag;



typedef enum cfiMeshType
{
	unknownMeshType,
	uniformQuadCornerType,		/* vertices at corners of quads     */
	uniformQuadCenterType		/* vertices at centers of quads     */

} cfiMeshType;



typedef enum cfiMeshSamplingMethod
{
	subNearestNeighbor,		/* subsample nearest neighbor by 2  */
	nearestNeighbor,		/* use the nearest neighbor         */
	linearInterp			/* linearly interpolate neighbors   */

} cfiMeshSamplingMethod;



typedef struct cfiMesh
{
	cfiMeshType  meshType;		/* type of mesh geometry            */
	cfiMatrix   *meshCoords;	/* encodes the mesh vertex coords   */
	cfiArray    *meshValues;	/* data values at the mesh vertices */

	/*
	 * The following cfiMeshTypes are supported.
	 *
	 * uniformQuadCornerType/uniformQuadCenterType
	 * -------------------------------------------
	 * An arbitrarility-dimensioned quadrilateral mesh whose
	 * 2-D cross-sections are planar parallelograms, e.g.
	 *
	 *	+---------+---------+---------+
	 *	 \         \         \         \
	 *	  +---------+---------+---------+
	 *	   \         \         \         \
	 *	    +---------+---------+---------+.
	 *
	 * The mesh vertices can be located either at the corners
	 * or at the centers of the volume elements bounded by the
	 * parallelograms.  Each mesh axis can have its own vertex
	 * spacing.
	 */

} cfiMesh;



/************************************************************************
 *									*
 * prototypes for functions in cfiMesh.c				*
 *									*
 ************************************************************************/

extern cfiErrorStatus	mapUniQuadMeshValues(
				const cfiMesh *const,
				const cfiMeshSamplingMethod,
				cfiMesh *const,
				cfiMatrix *const
			);
extern cfiErrorStatus	getUniQuadMeshIndexTransform
			(
				const cfiMesh *const,
				const cfiMesh *const,
				cfiMatrix *const
			);
extern void		sampleMeshValue(
				const cfiMeshSamplingMethod,
				const cfiCounter,
				const cfiCounter *const,
				const cfiScalar *const,
				const cfiDataType,
				void *,
				cfiCounter,
				const cfiCounter *,
				cfiScalar *const,
				cfiScalar *const
			);
extern cfiMesh *	allocateMesh(
				const cfiMeshType,
				const cfiCounter,
				const cfiDataType,
				const cfiCounter,
				const cfiCounter *const
			);
extern void		deallocateMesh(
				cfiMesh **const
			);
extern cfiErrorStatus	writeMesh(
				const cfiMesh *const,
				FILE *const
			);
extern cfiMeshStatus	checkMesh(
				const cfiMesh *const
			);
extern cfiErrorStatus	getUniQuadMeshOrigin(
				const cfiMesh *const,
				cfiScalar *
			);
extern cfiErrorStatus	setUniQuadMeshOrigin(
				const cfiScalar *,
				cfiMesh *const
			);
extern cfiErrorStatus	getUniQuadMeshAxis(
				const cfiMesh *const,
				const cfiDataAxis,
				cfiScalar *
			);
extern cfiErrorStatus	setUniQuadMeshAxis(
				const cfiScalar *,
				cfiMesh *const,
				const cfiDataAxis
			);
extern cfiMeshType	getMeshType(
				const cfiMesh *const
			);
extern cfiCounter	getMeshCoordDimension(
				const cfiMesh *const
			);
extern cfiCounter	getMeshAxisCount(
				const cfiMesh *const
			);
extern cfiArray *	getMeshDataArrayPtr(
				const cfiMesh *const
			);



#ifdef __cplusplus
}
#endif

#endif /* cfiMesh_HEADER */
