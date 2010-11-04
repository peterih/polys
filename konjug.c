/*
 * This module contains the source for a conjugate gradients algorithm.
 * Written by Soren Balling Engelsen, INRA-93.
 */
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include "polys.h"
#include "extern.h"

/* */
double
funk(double x[], double df[])
{
   register int  i;
   double        Vpot;

   for (i=0; i<M.nat; i++) {
      A[i].pos.x = x[3*i+1];
      A[i].pos.y = x[3*i+2];
      A[i].pos.z = x[3*i+3];
   }
   Vpot = potener(0);
   for (i=0; i<M.nat; i++) {
      df[3*i+1] = -F[i].x;
      df[3*i+2] = -F[i].y;
      df[3*i+3] = -F[i].z;
   }
   
   return Vpot;
}


#define FREEALL free_dvector(fm,1,n);free_dvector(h,1,n);free_dvector(g,1,n);
/* */
void
konjug(BOOLEAN silent, int n, double p[], int itmax, double ftol, 
       int *iter, double *fret, double (*funk)())
{
   register int  j;
   int      its, flag;
   double   gg, gam, fp, dgg, eps, la, lh;
   double   *g, *h, *fm;

   /* INITIALISATION */
   g  = dvector(1, n);
   h = dvector(1, n);
   fm = dvector(1, n);
   fp = (*funk)(p, fm);
   its = 0;
   flag = 0;
   gg = 0.0;
   lh = 0.0;
   for (j=1; j<=n; j++) 
   {  g[j] = -fm[j];
      h[j] = g[j];
      lh  += h[j]*h[j];
   }
   /* normalize search direction */
   lh = sqrt(lh);
   la = TOL;
   *fret = fp;
   if (!silent)
      printf("Iteration no. %d [%2d]  Epot = %.7lf kcal/mol  dS = %.7lf\n", 
              its, flag, *fret, lh*la);

   /* ITERATION LOOP */
   for (its=1; its<=itmax; its++) 
   {  *iter = its;
      eps = RODEPS;

      /* LINE SEARCH */
      linsch(n, &flag, &eps, &la, fret, fm, h, p, funk);
      if (!silent)
         printf("Iteration no. %d [%2d]  Epot = %.7lf kcal/mol  dS = %.7lf\n", 
                 its, flag, *fret, lh*la);

      /* UPDATE COORDINATE VECTOR AND FIND CONJUGATE DIRECTION */
      dgg = gg = 0.0;
      for (j=1; j<=n; j++) 
      {   p[j] += la*h[j];
          gg   += g[j]*g[j];
/*        dgg  += fm[j]*fm[j];                 Fletcher-Reeves  */
          dgg  += (fm[j] + g[j])*fm[j];     /* Polak-Ribiere    */
      } 

      switch(flag)
      {  case -5: printf("WARNING: Linesearch interpolation failed\n");
                  FREEALL;
                  return;
         case -4: printf("WARNING: Rounding errors ruined the interpolation\n");
                  break;
         case -3: printf("WARNING: Iterands very close in interpooler\n");
                  break;
         case -2: printf("WARNING: Linesearch extrapolation failed\n");
                  FREEALL;
                  return;
         case -1: printf("WARNING: Direction of line search is uphill\n");
                  FREEALL;
                  return;
         case  0: break;
      }

      if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPSILON)) 
      {  FREEALL
         printf("\nCONJUG: Exiting on function value criteria\n");
         return;
      }
      fp = *fret;
 
      if ((gg/n) < TOL) 
      /* the gradient  */
      {  FREEALL
         printf("\nCONJUG: Exiting succesfully on gradient criterium\n");
         return;
      }

      /* update search direction */
      lh = 0.0;
      gam = dgg/gg;
      for (j=1; j<=n; j++) 
      {  g[j] = -fm[j];
         h[j] = g[j]+gam*h[j];
         lh  += h[j]*h[j];
      }
      lh = sqrt(lh);

   }
   printf("\nCONJUG: Maximum number of iterations exceeded\n");
}
#undef FREEALL

/* */
void cartemini(BOOLEAN silent, int itmax, double ftol)
{
   register int  i;
   int      ndim, iter;
   double   fret, *p, rmsg;

   ndim = 3*M.nat;
   p = dvector(1, ndim);
   for (i=0; i<M.nat; i++)
   {  p[3*i+1] = A[i].pos.x;
      p[3*i+2] = A[i].pos.y;
      p[3*i+3] = A[i].pos.z;
   }
   rmsg = 0.0;
   nfun = 0;
   konjug(silent, ndim, p, itmax, ftol, &iter, &fret, funk);
   for (i=0; i<M.nat; i++)
   {  A[i].pos.x = p[3*i+1];
      A[i].pos.y = p[3*i+2];
      A[i].pos.z = p[3*i+3];
      rmsg += F[i].x*F[i].x;
      rmsg += F[i].y*F[i].y;
      rmsg += F[i].z*F[i].z;
   }
   rmsg = sqrt(rmsg/(3.0*M.nat));
   printf("CONJUG: Number of iterations: %10d\n", iter);
   printf("CONJUG: Number of function calls: %6d\n", nfun);
   printf("CONJUG: Func. value at solution: %15.7lf\n", fret);
   printf("CONJUG: RMS-gradient at solution: %14.7lf\n\n", rmsg);
   free_dvector(p, 1, ndim);

}
