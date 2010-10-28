/**********************************************************************
  "Vcon.c"

  This file contains the subroutines for calculation of the harmonic 
torsional constrain potential as well as its derivatives vrt. 
cartesians (F).

  Written by Soren Balling Engelsen, INRA-93.
**********************************************************************/
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "polys.h"
#include "extern.h"


/**********************************************************************
    Subroutine which will read in parameters and place them on top of
the torsional angle array (P).
The phase parameter is used as the constrain value.
The recommended constrain force constant is 500 kcal/mol.
         Ai - Aj - Ak - Al 
**********************************************************************/
void read_cp(int nphi, int *ncon)

{  FILE   *fp;
   int    i, torno, nocon;
   double Kcon, phicon, phi;
   char   line[MAXLEN];

   nocon = 0;
   if ((fp = fileopen("CONSTRAIN", "rt")) != NULL)
   {  get_line(fp, line, MAXLEN);
      printf("CONSTRAIN POTENTIAL -> %s\n", line);
      get_line(fp, line, MAXLEN);
      sscanf(line, "%d", &nocon);
      putchar('\n');

      for (i=nphi; i<nphi+nocon; i++)
      {  get_line(fp, line, MAXLEN);
         sscanf(line, "%d%lf%lf", &torno, &Kcon, &phicon);
         P[i] = P[torno]; 
         P[i].K1    = Kcon; 
         P[i].fold  = 0;
         P[i].p1 = phicon*DEGtoRAD;

         /*    -PI < cPHI < PI    */
         if (P[i].p1 > M_PI)
            P[i].p1 -= 2.0*M_PI;
         if (P[i].p1 < -M_PI)
            P[i].p1 += 2.0*M_PI;

         phi = torval(A[P[torno].from].pos, A[P[torno].over1].pos, 
                      A[P[torno].over2].pos, A[P[torno].to].pos);

         printf("Torsional Constraint no. %d :\n", i-nphi+1);
         printf("\tno. %3d  [%6.1lf kcal/mol]  PHI = %6.1lf   phi = %6.1lf\n",
                       torno, Kcon, phicon, phi);
      }
      fileclose("CONSTRAIN", fp);
   }
   printf("%d torsional constraints have been included\n\n", nocon);
   *ncon = nocon;

} /* End of read_cp */
 

/**********************************************************************
**********************************************************************/
double Vcon(int nphi, int ncon)

{  register int i;  
   BOOLEAN  deriv;
   int      ai, aj, ak, al;
   double   phi, phicon, phidif, pdsq, Kcon, w, V, dV;
   double   lT1, lT2, lTx;
   Vector3  RIJ, RJK, RKL, T1, T2, Tx;
   Vector3  dT1, dT2, dRIJ, dRJK, dRKL;

   V = 0.0;
   for (i=nphi; i<nphi+ncon; i++)
   {  /* get atom numbers */
      ai = P[i].from;
      aj = P[i].over1;
      ak = P[i].over2;
      al = P[i].to;
      phicon = P[i].p1;
      Kcon   = P[i].K1;

      /* get difference vectors */
      RIJ.x = A[ai].pos.x - A[aj].pos.x;
      RIJ.y = A[ai].pos.y - A[aj].pos.y;
      RIJ.z = A[ai].pos.z - A[aj].pos.z;

      RJK.x = A[aj].pos.x - A[ak].pos.x;
      RJK.y = A[aj].pos.y - A[ak].pos.y;
      RJK.z = A[aj].pos.z - A[ak].pos.z;

      RKL.x = A[al].pos.x - A[ak].pos.x;
      RKL.y = A[al].pos.y - A[ak].pos.y;
      RKL.z = A[al].pos.z - A[ak].pos.z;

      /* get crossproduct vectors */
      crossprod(RIJ, RJK, &T1);
      crossprod(RKL, RJK, &T2);

      /* normalize crossproduct vectors */
      lT1 = normalize(&T1);
      lT2 = normalize(&T2);
      if (lT1*lT2 < RODTOL) 
      {  printf("Torsion number %d is almost linear\n", i);
         w = 0.0;
         phi = 0.0;
         lTx = 0.0;
         deriv = FALSE;
      }
      else
      {  /* get the cosine relation */
         w = dotprod(T1, T2);
         if (w > 1.0)
            w = 1.0;
         if (w < -1.0)
            w = -1.0;
         phi = acos(w);
         deriv = TRUE;
      }
       
      /* determination of sign for phi */
      crossprod(T1, T2, &Tx);
      lTx = dotprod(Tx, Tx);
      lTx = sqrt(lTx);
      if (dotprod(Tx, RJK) > 0.0)
      {  phi = -phi;
         lTx = -lTx;
      }

      /*********************************************
      **    V = Kcon * (phi - phicon)**2          **
      *********************************************/

      phidif = (phi - phicon);
      if (phidif > M_PI)
         phidif -= 2.0*M_PI;
      if (phidif < -M_PI)
         phidif += 2.0*M_PI;
      
      pdsq  = phidif * phidif;
      if ((fabs(phicon) < TOL) AND (fabs(phidif) < 0.1))
      {  V += Kcon * pdsq;
         dV = -2.0 * Kcon * (1.0 + pdsq/6.0);
      }
      else 
      if (((M_PI - fabs(phicon)) < TOL) AND (fabs(phidif) < 0.1))
      {  V += Kcon * pdsq;
         dV = 2.0 * Kcon * (1.0 + pdsq/6.0);
      }
      else
      {  if (fabs(lTx) < RODEPS)
         {   if (lTx < 0.0)
               lTx = -RODEPS;
             else
               lTx = RODEPS;
             printf("lTx   : %.10lf\n", lTx);
          }
          V += Kcon * pdsq;
          dV = -2.0 * Kcon * phidif / lTx;
      }
         
      if (deriv)
      {  /* calculate dV/cosx */       

         dT1.x = (T2.x - w*T1.x)/lT1;
         dT1.y = (T2.y - w*T1.y)/lT1;
         dT1.z = (T2.z - w*T1.z)/lT1;

         dT2.x = (T1.x - w*T2.x)/lT2;
         dT2.y = (T1.y - w*T2.y)/lT2;
         dT2.z = (T1.z - w*T2.z)/lT2;

         crossprod(RJK, dT1, &dRIJ);
         crossprod(RJK, dT2, &dRKL);

         dRJK.x = dT1.y*RIJ.z - dT1.z*RIJ.y + dT2.y*RKL.z - dT2.z*RKL.y;
         dRJK.y = dT1.z*RIJ.x - dT1.x*RIJ.z + dT2.z*RKL.x - dT2.x*RKL.z;
         dRJK.z = dT1.x*RIJ.y - dT1.y*RIJ.x + dT2.x*RKL.y - dT2.y*RKL.x;

         /* calculate forces  */
         F[ai].x -= dRIJ.x * dV;
         F[ai].y -= dRIJ.y * dV;
         F[ai].z -= dRIJ.z * dV;

         F[aj].x += (dRIJ.x - dRJK.x) * dV;
         F[aj].y += (dRIJ.y - dRJK.y) * dV;
         F[aj].z += (dRIJ.z - dRJK.z) * dV;

         F[ak].x += (dRJK.x + dRKL.x) * dV;
         F[ak].y += (dRJK.y + dRKL.y) * dV;
         F[ak].z += (dRJK.z + dRKL.z) * dV;

         F[al].x -= dRKL.x * dV;
         F[al].y -= dRKL.y * dV;
         F[al].z -= dRKL.z * dV;
      }
   }

   return (V);

}  /* End of Vcons */


/* End of file */
