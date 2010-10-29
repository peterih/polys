/**********************************************************************
   "Init.c"

   This module contains the subroutine "main" to program polys.

   Written by Soren Balling Engelsen, INRA-93.
**********************************************************************/
#include "defs.h"
#include "typs.h"
#include "extern.h"


/**********************************************************************
	Initializes global data variables.
**********************************************************************/
void Initialize(void)

{  int   i;
 
   /* timer variables */
   tbo = 0.0;
   tv = 0.0;
   tt = 0.0;
   trot = 0.0;
   tnb = 0.0;
   tbuild = 0.0;
   ttest = 0.0;
   tpara = 0.0;
   tint = 0.0;
   tint = 0.0;

   /* molecule variables */
   da = 0;
   M.nat   = 0;
   M.nbond = 0;
   M.ntheta= 0;
   M.nphi  = 0;
   M.noop  = 0;
   M.nres  = 1;
   M.nrot  = 0;
   M.Eb = 0.0;
   M.Et = 0.0;

   /* atom variables */
   for ( i=0; i<MAXATOM; i++)
   {  A[i].res   = 0;
      A[i].pos.x = 0.0;
      A[i].pos.y = 0.0;
      A[i].pos.z = 0.0;
      A[i].chg   = 0.0;
      A[i].ztyp  = 0;
      F[i].x = 0.0;
      F[i].y = 0.0;
      F[i].z = 0.0;
      DA[i] = 0;
   }
   for ( i=0; i<MAXRES; i++)
   {  Res[i].ringtype = 0;
      Res[i].status = 0;
      Res[i].from = 0;
      Res[i].to = 0;
      Res[i].conres = 0;
      Res[i].theta = 999.9;
      Res[i].phi = 999.9;
      Res[i].psi = 999.9;
      Res[i].omega = 999.9;
   }
      
   if (debug) 
   {  printf("INIT: DIMENSIONS OF DATASTRUCTURES\n");
      printf("INIT:\tMAXATOM  : %d\n", MAXATOM);
      printf("INIT:\tMAXRES   : %d\n", MAXRES);
      printf("INIT:\tMAXTYPE  : %d\n", MAXTYPE);
      printf("INIT:\tMAXBOND  : %d\n", MAXBOND);
      printf("INIT:\tMAXTHETA : %d\n", MAXTHETA);
      printf("INIT:\tMAXPHI   : %d\n", MAXPHI);
      printf("INIT:\tMAXROT   : %d\n", MAXROT);
      printf("INIT:\tMAXOOP   : %d\n", MAXOOP);
      printf("INIT:\tEPSILON  : %.20f\n", EPSILON);
      printf("INIT:\tRODEPS   : %.20f\n", RODEPS);
      printf("INIT:\tTOL      : %.20f\n", TOL);
      printf("INIT:\tRODTOL   : %.20f\n\n", RODTOL);
      printf("INIT: SIZE OF DATASTRUCTURES\n");
      printf("INIT:\tchar  : %ld\n", sizeof(char));
      printf("INIT:\tshort : %ld\n", sizeof(short));
      printf("INIT:\tint   : %ld\n", sizeof(int));
      printf("INIT:\tlong  : %ld\n", sizeof(long));
      printf("INIT:\tfloat : %ld\n", sizeof(float));
      printf("INIT:\tdouble: %ld\n", sizeof(double));
      printf("INIT:\tA[1]  : %ld\n", sizeof(A[1]));
      printf("INIT:\tA[]   : %ld\n", sizeof(A));
      printf("INIT:\tRes[1]: %ld\n", sizeof(Res[1]));
      printf("INIT:\tRes[] : %ld\n", sizeof(Res));
      printf("INIT:\tB[1]  : %ld\n", sizeof(B[1]));
      printf("INIT:\tB[]   : %ld\n", sizeof(B));
      printf("INIT:\tT[1]  : %ld\n", sizeof(T[1]));
      printf("INIT:\tT[]   : %ld\n", sizeof(T));
      printf("INIT:\tP[1]  : %ld\n", sizeof(P[1]));
      printf("INIT:\tP[]   : %ld\n", sizeof(P));
      printf("INIT:\tRT[1] : %ld\n", sizeof(RT[1]));
      printf("INIT:\tRT[]  : %ld\n", sizeof(RT));
      printf("INIT:\tO[1]  : %ld\n", sizeof(O[1]));
      printf("INIT:\tO[]   : %ld\n", sizeof(O));
      printf("INIT:\tM     : %ld\n\n", sizeof(M));
   }
} /* End of initialize */


/* End of file */
