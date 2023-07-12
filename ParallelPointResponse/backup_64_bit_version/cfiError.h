/*
 * SCCS:  @(#)cfiError.h  1.6  6/7/95  11:33:52
 */

#ifndef cfiError_HEADER
#define cfiError_HEADER

#ifdef __cplusplus
extern "C" {
#endif



#include <stdio.h>
#include <math.h>



typedef unsigned int cfiErrorStatus;

typedef enum cfiErrorFlag
{
	cfiSuccess = 0x00,
	cfiWarning = 0x01,
	cfiFailure = 0x02

} cfiErrorFlag;



typedef unsigned int cfiDebugStatus;

typedef enum cfiDebugFlag
{
	cfiNoDebug       = 0x00,	/* no debugging info                */
	cfiMinDebug      = 0x01,	/* minimum of eight debug levels    */
	cfiProgressDebug = 0x02,	/* high level program progress info */
	cfiInOutDebug    = 0x04,	/* user/file input/output info      */
	cfiMemoryDebug   = 0x08,	/* memory management info           */
	cfiTraceDebug    = 0x10,	/* low level function progress info */
	cfiArgDebug      = 0x20,	/* do function argument checking    */
	cfiReturnDebug   = 0x40,	/* suspicious function return info  */
	cfiMaxDebug      = 0x80		/* maximum of eight debug levels    */

} cfiDebugFlag;

extern cfiDebugStatus g_debugStatus;	/* global debug flags               */



/************************************************************************
 *									*
 * message macros							*
 *									*
 ************************************************************************/

#define MSG_LENGTH		(256)	/* for compound messages            */

#define SHOW_LINE_IN_MSGS	(  1)
#if ( SHOW_LINE_IN_MSGS )

#define WARNING_MSG(msg)						\
	(void)fprintf(stderr, "warning(%s,%d): %s\n",			\
		      __FILE__, __LINE__, msg)

#define FAILURE_MSG(msg)						\
	(void)fprintf(stderr, "FAILURE(%s,%d): %s\n",			\
		      __FILE__, __LINE__, msg)

#define DEBUG_MSG(flag, msg)						\
	if ( flag & g_debugStatus )					\
		(void)fprintf(stderr, "debug%d(%s,%d): %s\n",		\
			      (int)((log(2 * flag) / log(2)) + 0.5),	\
			      __FILE__, __LINE__, msg)

#else

#define WARNING_MSG(msg)						\
	(void)fprintf(stderr, "warning(%s): %s\n",			\
		      __FILE__, msg)

#define FAILURE_MSG(msg)						\
	(void)fprintf(stderr, "FAILURE(%s): %s\n",			\
		      __FILE__, msg)

#define DEBUG_MSG(flag, msg)						\
	if ( flag & g_debugStatus )					\
		(void)fprintf(stderr, "debug%d(%s): %s\n",		\
			      (int)((log(2 * flag) / log(2)) + 0.5),	\
			      __FILE__, msg)

#endif



#ifdef __cplusplus
}
#endif

#endif /* cfiError_HEADER */
