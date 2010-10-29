/*
 * This file contains the subroutines for calculation of torsional angle
 * potential energy as well as its derivatives vrt. cartesians (F).
 * Written by Soren Balling Engelsen, INRA-93.
 */
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "polys.h"
#include "extern.h"

static int peftype;

/*
 * Function which compares to torsion types to see weather they fits.
 * Ai - Aj - Ak  - Al
 */
BOOLEAN
torcomp(int it, int jt, int kt, int lt, Tors Pi)
{
   int   tf, to1, to2, tt;

   tf  = A[Pi.from].ztyp;
   to1 = A[Pi.over1].ztyp;
   to2 = A[Pi.over2].ztyp;
   tt  = A[Pi.to].ztyp;
   if ((tf == it) OR (it == 9999))
     if (to1 == jt)
       if (to2 == kt)
         if ((tt == lt) OR (lt == 9999))
           return (TRUE);
   if ((tt == it) OR (it == 9999))
     if (to2 == jt)
       if (to1 == kt)
         if ((tf == lt) OR (lt == 9999))
           return (TRUE);
   return (FALSE);
}

/*
 * Subroutine which will read in parameters and place them in the 
 * torsional angle array (P). The routine avoids an extra parameter array.
 * Ai - Aj - Ak - Al 
 */
void
read_tp(int nphi)
{
   FILE   *fp;
   int    i, j, nopara, it, jt, kt, lt;
   int    fold, count = 0;
   double K1, p1, K2, p2, K3, p3;
   char   line[MAXLEN];
   char   *pardir = getenv("PARAM_DIR");

   strcpy(line,pardir);
   strcat(line,"/TPARAM.PF");
   if ((fp = fileopen(line, "rt")) != NULL)
   {  get_line(fp, line, MAXLEN);
      printf("TORSION POTENTIAL -> %s\n", line);
      get_line(fp, line, MAXLEN);
      sscanf(line, "%d%d", &nopara, &peftype);

      switch(peftype)
      {  case 0: printf("POTENTIAL IGNORED");
                 break;
         case 1: printf("SIMPLE PITZER POTENTIAL");
                 break;
         case 2: printf("PITZER POTENTIAL WITH PHASE");
                 break;
         case 3: printf("FOURIER SERIE WITH PHASE");
                 break;
         default: printf("POTENTIAL IGNORED");
      }
      putchar('\n');
      putchar('\n');

      for (i=0; i<nopara; i++)
      {  get_line(fp, line, MAXLEN);
         switch(peftype)
         { case 1: /* Simple pitzer */
                   sscanf(line, "%d%d%d%d%lf%d", 
                                &it, &jt, &kt, &lt, &K1, &fold);
                             p1 = 0.0;
                   K2 = 0.0; p2 = 0.0;
                   K3 = 0.0; p3 = 0.0;
                   break;
           case 2: /* Pitzer with phase */
                   sscanf(line, "%d%d%d%d%lf%d%lf", 
                                &it, &jt, &kt, &lt, &K1, &fold, &p1);
                   K2 = 0.0; p2 = 0.0;
                   K3 = 0.0; p3 = 0.0;
                   break;
           case 3: /* Fourier with phase */
                   sscanf(line, "%d%d%d%d%lf%lf%lf%lf%lf%lf", 
                     &it, &jt, &kt, &lt, &K1, &p1, &K2, &p2, &K3, &p3);
                   fold = 0;
                   break;
         }
         for (j=0; j< nphi; j++)
         {  if (torcomp(it, jt, kt, lt, P[j]))
            { P[j].fold = fold;
              P[j].K1 = K1/2.0;
              if (fabs(p1) < 0.01)
                 P[j].p1 = 1.0;
              else
                 P[j].p1 = -1.0;
              P[j].K2 = K2/2.0;
              if (fabs(p2) < 0.01)
                 P[j].p2 = 1.0;
              else
                 P[j].p2 = -1.0;
              P[j].K3 = K3/2.0;
              if (fabs(p3) < 0.01)
                 P[j].p3 = 1.0;
              else
                 P[j].p3 = -1.0;
              count++;
            }
         }
      }
      putchar('\n');
      fileclose("TPARAM.PF", fp);
   }
   if (count < nphi)
   {  printf("WARNING: MISSING TORSION PARAMETERS\n");
      printf("WARNING: Missing parameters to %d torsional angles\n",
                       nphi - count);
      for (j=0; j< nphi; j++)
         if (P[j].fold == 0)
            printf("%d: %d - %d - %d - %d\n", j, A[P[j].from].ztyp, 
                   A[P[j].over1].ztyp, A[P[j].over2].ztyp, A[P[j].to].ztyp);
      printf("WARNING: MISSING TORSION PARAMETERS\n\n");
   }
} /* End of read_tp */
 

/* */
double
Vphi(int nphi)
{
   register int i;  
   int      ai, aj, ak, al;
   double   tstart;
   double   K1, p1, K2, p2, K3, p3;
   double   phi, w, sqw, V, dV;
   double   lT1, lT2, lTx;
   Vector3  RIJ, RJK, RLK, T1, T2, Tx;
   Vector3  dT1, dT2, dRIJ, dRJK, dRLK;
   BOOLEAN  deriv;

   tstart = cput();
   V = 0.0;
   for (i=0; i<nphi; i++)
   {  /* get atom numbers */
      ai = P[i].from;
      aj = P[i].over1;
      ak = P[i].over2;
      al = P[i].to;

      /* get difference vectors */
      RIJ.x = A[ai].pos.x - A[aj].pos.x;
      RIJ.y = A[ai].pos.y - A[aj].pos.y;
      RIJ.z = A[ai].pos.z - A[aj].pos.z;

      RJK.x = A[aj].pos.x - A[ak].pos.x;
      RJK.y = A[aj].pos.y - A[ak].pos.y;
      RJK.z = A[aj].pos.z - A[ak].pos.z;

      RLK.x = A[al].pos.x - A[ak].pos.x;
      RLK.y = A[al].pos.y - A[ak].pos.y;
      RLK.z = A[al].pos.z - A[ak].pos.z;

      /* get crossproduct vectors */
      crossprod(RIJ, RJK, &T1);
      crossprod(RLK, RJK, &T2);

      /* normalize crossproduct vectors */
      lT1 = normalize(&T1);
      lT2 = normalize(&T2);

      if (lT1*lT2 < 0.01) 
      {  printf("Torsion number %d is almost linear\n", i);
         deriv = FALSE;
      }
      else
         deriv = TRUE;

      /* get the cosine relation */
      w = dotprod(T1, T2);
      if (w > 1.0)
         w = 1.0;
      if (w < -1.0)
         w = -1.0;
      phi = acos(w);
 
      /* determination of sign for phi */
      crossprod(T1, T2, &Tx);
      lTx = length(Tx);
      if (dotprod(Tx, RJK) > 0.0)
      {  phi = -phi;
         lTx = -lTx;
      }

      sqw = w*w;
      K1 = P[i].K1; p1 = P[i].p1;
      K2 = P[i].K2; p2 = P[i].p2;
      K3 = P[i].K3; p3 = P[i].p3;
      switch(P[i].fold)
      {  case  0: 
         /**************************************************************/
         /*   V =   K1/2*[1 + p1*cos(a)]                               */
         /*       + K2/2*[1 + p2*cos(2a)]                              */
         /*       + K3/2*[1 + p3*cos(3a)]                              */
         /**************************************************************/
                  V  += K1*(1.0 + p1*w);
                  V  += K2*(1.0 + p2*2.0*sqw - p2);   
                  V  += K3*(1.0 + p3*4.0*sqw*w - p3*w*3.0);
                  dV  = p1*K1;
                  dV += p2*K2*4.0*w;
                  dV += p3*K3*(12.0*sqw - 3.0);
                  break;
         case  3: 
         /**************************************************************/
         /*   cos(3a) = 4*cos(a)**3 - 3*cos(a)                         */
         /**************************************************************/
                  V  += K1*(1.0 + 4.0*sqw*w - 3.0*w);
                  dV  = K1*(12.0*sqw - 3.0);
                  break;
         case -2: 
         /**************************************************************/
         /*   cos(2a) = 2*cos(a)**2 - 1                                */
         /**************************************************************/
                  V  += K1*2.0*(1.0 - sqw);   
                  dV  = -4.0*K1*w;
                  break;
         case -4: 
         /**************************************************************/
         /* fourth order polynomia  - EXPERIMENTAL                     */
         /**************************************************************/
                  V  += K1*sqw*sqw;   
                  dV  = 4.0*K1*sqw*w;
                  break;
         case  4: 
         /**************************************************************/
         /*   cos(4a) = 8*cos(a)**4 - 8*cos(a)**2 + 1                  */
         /**************************************************************/
                  V  += K1*(2.0 + 8.0*sqw*sqw - 8.0*sqw);
                  dV  = K1*8.0*(2.0*sqw*w - w);
                  break;
         case  6: 
         /**************************************************************/
         /*   cos(6a) = 32*cos(a)**6 - 48*cos(a)**4 + 18*cos(a)**2 - 1 */
         /**************************************************************/
                  V  += K1*2.0*(16.0*sqw*sqw*sqw - 24.0*sqw*sqw + 9.0*sqw);
                  dV  = K1*2.0*(96.0*sqw*sqw*w - 96.0*sqw*w + 18.0*w);
                  break;
         default: 
                  printf("Special %d-fold torsional potential not supported\n",
                          P[i].fold);
                  V  += 0.0;
                  dV  = 0.0;
      }

      if (deriv)
      {  /* calculate dV/cosx */       
         dT1.x = (T2.x - w*T1.x)/lT1;
         dT1.y = (T2.y - w*T1.y)/lT1;
         dT1.z = (T2.z - w*T1.z)/lT1;

         dT2.x = (T1.x - w*T2.x)/lT2;
         dT2.y = (T1.y - w*T2.y)/lT2;
         dT2.z = (T1.z - w*T2.z)/lT2;

         dRIJ.x = dT1.z*RJK.y - dT1.y*RJK.z;
         dRIJ.y = dT1.x*RJK.z - dT1.z*RJK.x;
         dRIJ.z = dT1.y*RJK.x - dT1.x*RJK.y;

         dRJK.x = dT1.y*RIJ.z - dT1.z*RIJ.y + dT2.y*RLK.z - dT2.z*RLK.y;
         dRJK.y = dT1.z*RIJ.x - dT1.x*RIJ.z + dT2.z*RLK.x - dT2.x*RLK.z;
         dRJK.z = dT1.x*RIJ.y - dT1.y*RIJ.x + dT2.x*RLK.y - dT2.y*RLK.x;

         dRLK.x = dT2.z*RJK.y - dT2.y*RJK.z;
         dRLK.y = dT2.x*RJK.z - dT2.z*RJK.x;
         dRLK.z = dT2.y*RJK.x - dT2.x*RJK.y;

         /* calculate forces  */
         F[ai].x -= dRIJ.x * dV;
         F[ai].y -= dRIJ.y * dV;
         F[ai].z -= dRIJ.z * dV;

         F[aj].x += (dRIJ.x - dRJK.x) * dV;
         F[aj].y += (dRIJ.y - dRJK.y) * dV;
         F[aj].z += (dRIJ.z - dRJK.z) * dV;

         F[ak].x += (dRJK.x + dRLK.x) * dV;
         F[ak].y += (dRJK.y + dRLK.y) * dV;
         F[ak].z += (dRJK.z + dRLK.z) * dV;

         F[al].x -= dRLK.x * dV;
         F[al].y -= dRLK.y * dV;
         F[al].z -= dRLK.z * dV;
      }
   }

   tt += cput() - tstart;
   return (V);
}  /* End of Vphi */
