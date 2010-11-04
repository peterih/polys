/*
 * This module contains the code for a strict line search algorithm,
 * using cubic interpolation.
 * Written by Soren Balling Engelsen, INRA-93.
 */
# include <math.h>
# include "defs.h"

#define  EXTRA  20
#define  INTER  20
/*
 * INPUT:
 * n       : the dimension of the function to minimize
 * *flag   :
 * *eps    :
 * *la     :
 * *y      :
 * fm[]    :
 * h[]     :
 * x[]     :
 * func    :
 *
 * OUTPUT:
 * *flag   :
 * *eps    :
 * *la     :
 * *y      :
 * fm[]    :
 * h[]     :
 * x[]     :
 * func    :
 */
void
linsch(int n, int *flag, double *eps, double *la, double *y, 
       double fm[], double h[], double x[], double (*func)())
{
   int    i, eit, iit;
   double t1, t2, t3, tg, f1, f2, f3, tm1, fm1, fm2, fm3, eta, gam, l;
   double *xa, *dvector();

   /* initialize */
   *flag = 0;      /* OK */
   xa = dvector(1,n);
   t2 = 0.0;
   l  = *la;
   f2 = *y;
   fm2 = fm[1]*h[1];
   for (i=2; i<=n; i++)
      fm2 += fm[i]*h[i];

   /* search direction "h" is uphill */
   if (fm2 > 0.0)
   {  *flag = -1; 
      *la = 0.0;
      *eps = 0.0;
      for (i=1; i<=n; i++)
         h[i] = 0.0;
      free_dvector(xa,1,n);
      return;
   }

   /* EXTRAPOLATION */
   eit = 1;
   do 
   {  t1 = t2;
      f1 = f2;
      fm1 = fm2;
      t2 += l;
      for (i=1; i<=n; i++)
         xa[i] = x[i] + t2*h[i];
      f2 = (*func)(xa, fm);
      eit++;
      fm2 = fm[1]*h[1]; 
      for (i=2; i<=n; i++)
         fm2 += fm[i]*h[i];
      l *= 1.62;
      if (fm2 > 0.0)          /* jump out when uphill */
         break;
   } while (eit < EXTRA);

   /* no finite solution */
   if (eit == EXTRA)
   {  *flag = -2;
      *la = t2;
      *y = f2;
      *eps = t2 - t1;
      free_dvector(xa,1,n);
      return;
   }

   /* INTERPOLATION */
   iit = 0;
   tg = -t2;
   do
   {  eta = (f1 - f2)*3.0/(t2 - t1) + fm1 + fm2;
      gam = sqrt(eta*eta - fm1*fm2);
      t3  = t1 + (t2 - t1)*(1.0 - (fm2 + gam - eta)/(fm2 - fm1 + gam*2.0));

      /* iterands very close */ 
      if (fabs(tg-t3) < (t2+t1)*EPSILON)
      {  *flag = -3;
         *la = t3;
         *y = f3;
         *eps = t2 - t1;
         free_dvector(xa,1,n);
         return;
      }
      
      if ((t3 < t1) OR (t3 > t2))
      {  *flag = -4;
         *la = t3;
         *y = f3;
         *eps = t2 - t1;
         free_dvector(xa,1,n);
         return;
      }

      /* renewal */
      for (i=1; i <=n; i++)
         xa[i] = x[i] + t3*h[i];
      f3 = (*func)(xa, fm);
      iit++;
      fm3 = fm[1]*h[1];
      for (i=2; i<=n; i++)
         fm3 += fm[i]*h[i]; 
      tg = t3;

      if ((f3 < f2) AND (f3 < f1))
         break;

      if (fm3 <= 0.0)
      /* right hand interval part */
      {  t1 = t3;
         f1 = f3;
         fm1 = fm3;
      }
      else
      /* left hand interval part */
      {  t2 = t3;
         f2 = f3;
         fm2 = fm3;
      }
   } while (iit<INTER);

   /* no interpolate */
   if (iit == INTER)
   {  *flag = -5;
      *la = t2;
      *y = f2;
      *eps = t2 - t1;
      free_dvector(xa,1,n);
      return;
   }

   *flag = eit+iit;
   *y = f3;
   *la = t3;
   *eps = t2 - t1;
   free_dvector(xa,1,n);
   return;
}
