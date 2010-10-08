/*****************************************************************************
	"defs.h" 
	This file contains a number of convinient definitions of global 
	constants and logical operators.
	SBE93-INRA
*****************************************************************************/

/* FIRST some general convinient definitions */

#ifndef  MAXLEN
#define  MAXLEN   80
#endif

#ifndef  NULL
#define  NULL     ((void *) 0)
#endif

#ifndef  M_PI
#define  M_PI     3.14159265358979323846
#endif

#define  RADtoDEG 180.0/M_PI
#define  DEGtoRAD M_PI/180.0

#define  GOAL     4000000L

#ifndef  MAX
#define  MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#ifndef  MIN
#define  MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#ifndef  INRANGE
#define  INRANGE(x,y,v) ((v) >= (x) && (v) <= (y))
#endif

#ifndef  SWAP
#define  SWAP(tmp,x,y)  (tmp = x, x = y, y = tmp)
#endif

#ifndef  BOOLEAN
#define  BOOLEAN    short
#define  BOOLVAL(x) (x) ? "TRUE" : "FALSE"
#define  OR       ||
#define  AND      &&
#define  NOT      !=
#define  FALSE    0
#define  TRUE     1
#endif

#define  STREQ(x,y) (strcmp(x,y) == 0)
#define  STRLT(x,y) (strcmp(x,y) < 0)
#define  STRGT(x,y) (strcmp(x,y) > 0)

#define  ALLOC(x)   ((x *) malloc((unsigned) sizeof(x)))

/* THEN Molecular Mechanics specific definitions */

#define  EPSILON  1.0E-16
#define  RODEPS   1.0E-8
#define  TOL      1.0E-6
#define  RODTOL   1.0E-3

#define  MAXTINT  10000
#define  MAXATOMR   200
#define  MAXCOL    1000

#ifdef __TURBOC__
#define  MAXTYPE   2000
#define  MAXRES      25
#define  MAXIMG     108
#define  MAXATOM    500
#define  MAXBOND    600
#define  MAXTHETA  1000
#define  MAXPHI    1500
#define  MAXROT     200
#define  MAXOOP      10
#else
#define  MAXTYPE   2000
#define  MAXRES    2500
#define  MAXIMG     108
#define  MAXATOM  50000
#define  MAXBOND  60000
#define  MAXTHETA 100000
#define  MAXPHI   150000
#define  MAXROT   20000
#define  MAXOOP     200
#endif
/* End of file */
