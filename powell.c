/**********************************************************************
  "powell.c"

  This module contains the nonderivative Powell rutine for minimizing 
a function. See Numerical Recipes p. 299.

  Adapted for POLYS by Soren Balling Engelsen, INRA-93.
**********************************************************************/
#include <stdio.h>
#include <math.h>
#include "polys.h"
#include "extern.h"

static double sqrarg;
#define SQR(a) (sqrarg=(a),sqrarg*sqrarg)
/**********************************************************************
   option is used for printing:

	-1: No printout
	 0: Energy printout (phimin).
	 1: dn2 printout (helix).
**********************************************************************/
void powell(BOOLEAN option, double *p, double **xi, int ndim,
            double ftol, int *iter, double *fret, double (*func)())

{  int    i, ibig, j;
   double t,fptt,fp,del,new,tole=0.02;
   double *pt,*ptt,*xit;

   /* Initializations */
   pt=dvector(1,ndim);
   ptt=dvector(1,ndim);
   xit=dvector(1,ndim);
   *fret=(*func)(p);
   fp=(*fret);
   for (j=1;j<=ndim;j++) 
      pt[j]=p[j];
   switch (option)
   {  case -1 : printf("Iteration no. %d [%2d]  Epot = %.7lf kcal/mol",
                       "  dE = %.7lf\n",0, nfun, (*fret), fabs(fp-(*fret)));
      case  0 : break;
      case  1 : printf("Iteration no. %d [%2d]  dN2 = %.7lf\n",0,nfun,(*fret));
                print = TRUE;
                new=(*func)(p);
                print = FALSE;
                break;
   }

   /* Main loop */
   for (*iter=1;;(*iter)++) 
   {  fp=(*fret);
      ibig=0;
      del=0.0;
      for (i=1;i<=ndim;i++) 
      {  for (j=1;j<=ndim;j++) 
            xit[j]=xi[j][i];
         fptt=(*fret);
         linmin(p,xit,ndim,tole,fret,func);
/*         printf("%12.7f\n",(*fret));     */
         if (fabs(fptt-(*fret)) > del) 
         {  del=fabs(fptt-(*fret));
            ibig=i;
         }
      }
      switch (option)
      {  case -1 : printf("Iteration no. %d [%2d]  Epot = %.7lf kcal/mol",
                      "  dE = %.7lf\n",(*iter),nfun,(*fret),fabs(fp-(*fret)));
         case  0 : break;
         case  1 : printf("Iteration no. %d [%2d]  dN2 = %.7lf\n",(*iter),
                           nfun,(*fret));
                   print = TRUE;
                   new=(*func)(p);
                   print = FALSE;
                   break;
      }
      if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) 
      {  free_dvector(xit,1,ndim);
         free_dvector(ptt,1,ndim);
         free_dvector(pt,1,ndim);
         printf("\n\nPOWELL: Optimization succesfully accomplished\n");
         return;
      }
      if (*iter > maxit) 
      {  free_dvector(xit,1,ndim);
         free_dvector(ptt,1,ndim);
         free_dvector(pt,1,ndim);
         printf("\n\nPOWELL: Maximum number of iterations exceeded!!\n");
         return;
      }
      for (j=1;j<=ndim;j++) 
      {  ptt[j]=2.0*p[j]-pt[j];
         xit[j]=p[j]-pt[j];
         pt[j]=p[j];
      }
      fptt=(*func)(ptt);
      if (fptt < fp) 
      {  t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
         if (t < 0.0) 
         {  linmin(p,xit,ndim,tole,fret,func);
            for (j=1;j<=ndim;j++) 
               xi[j][ibig]=xit[j];
         }
      }
   }
} /* End of powell */
#undef SQR

/* End of file */
