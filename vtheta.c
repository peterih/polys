/**********************************************************************
  "Vtheta.c"

  This file contains the subroutines for calculation of valence angle
potential energy as well as its derivatives vrt. cartesians (F).

  Written by Soren Balling Engelsen, INRA-93.
**********************************************************************/
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "polys.h"
#include "extern.h"

static int peftype;

/**********************************************************************
    Function which creates a unique code for each type of valence angle.
         Ai - Aj - Ak 
**********************************************************************/
double anglecode(Angl Ti)

{  int  it, jt, kt;

   it = A[Ti.from].ztyp;
   jt = A[Ti.over].ztyp;
   kt = A[Ti.to].ztyp;
   if (kt > it)
     return (kt + jt * MAXTINT + it * MAXTINT * MAXTINT);
   else
     return (it + jt * MAXTINT + kt * MAXTINT * MAXTINT); 

} /* End of anglecode */

/**********************************************************************
    Function which creates a unique code for each type of valence angle.
         Ai - Aj - Ak 
**********************************************************************/
double acode(int it, int jt, int kt)

{
   if (kt > it)
     return (kt + jt * MAXTINT + it * MAXTINT * MAXTINT);
   else
     return (it + jt * MAXTINT + kt * MAXTINT * MAXTINT); 

} /* End of acode */

/**********************************************************************
    Subroutine which will read in parameters and place then in the 
valence angle array (T). The routine avoids and extra parameter array.
         Ai - Aj - Ak 
**********************************************************************/
void read_ap(int ntheta)

{  FILE   *fp;
   int    i, j, nopara, it, jt, kt;
   int    count = 0;
   double aicode, ajcode, Ktheta, t0;
   char   line[MAXLEN];
   char   *pardir = getenv("PARAM_DIR");

   strcpy(line,pardir);
   strcat(line,"/APARAM.PF");
   if ((fp = fileopen(line, "rt")) != NULL)
   {  getline(fp, line, MAXLEN);
      printf("ANGLE POTENTIAL -> %s\n", line);
      getline(fp, line, MAXLEN);
      sscanf(line, "%d%d", &nopara, &peftype);

      switch(peftype)
      {  case 0: printf("POTENTIAL IGNORED");
                 break;
         case 1: printf("HARMONIC POTENTIAL");
                 break;
         case 2: printf("ROVAP POTENTIAL");
                 break;
        default: printf("POTENTIAL IGNORED");
      }
      putchar('\n');
      putchar('\n');

      for (i=0; i<nopara; i++)
      {  getline(fp, line, MAXLEN);
         sscanf(line, "%d%d%d%lf%lf", &it, &jt, &kt, &Ktheta, &t0);
         aicode = acode(it, jt, kt);
         for (j=0; j< ntheta; j++)
         {  ajcode = anglecode(T[j]);
            if (ajcode == aicode)
            { T[j].Kt = Ktheta;
              T[j].t0 = t0*DEGtoRAD;
              count++;
            }
         }
      }
      putchar('\n');
      fileclose("APARAM.PF", fp);
   }
   if (count !=  ntheta)
   {  printf("WARNING: MISSING ANGLE PARAMETERS\n");
      printf("WARNING: Missing parameters to %d valence angles\n",
                       ntheta - count);
      for (j=0; j< ntheta; j++)
         if (T[j].t0 < 0.1)
            printf("%d: %d - %d - %d\n", j, A[T[j].from].ztyp,
                   A[T[j].over].ztyp, A[T[j].to].ztyp);
      printf("WARNING: MISSING ANGLE PARAMETERS\n\n");
   }

} /* End of read_ap */
 
/**********************************************************************
**********************************************************************/
double Vtheta(int ntheta)

{  register int i;  
   int      ai, aj, ak;
   double   tstart;
   double   theta, w, RIJsq, RKJsq, lRIJ, lRKJ, tt, V, dV, dtdw;
   Vector3  RIJ, RKJ, dtdRIJ, dtdRKJ;

   tstart=cput();
   V = 0.0;
   if (peftype != 0)
   for (i=0; i<ntheta; i++)
   {  /* get atom numbers */
      ai = T[i].from;
      aj = T[i].over;
      ak = T[i].to;
      /* get difference vectors */
      RIJ.x = A[ai].pos.x - A[aj].pos.x;
      RIJ.y = A[ai].pos.y - A[aj].pos.y;
      RIJ.z = A[ai].pos.z - A[aj].pos.z;
      RKJ.x = A[ak].pos.x - A[aj].pos.x;
      RKJ.y = A[ak].pos.y - A[aj].pos.y;
      RKJ.z = A[ak].pos.z - A[aj].pos.z;
      /* normalize difference vectors */
      RIJsq = dotprod(RIJ, RIJ);
      RKJsq = dotprod(RKJ, RKJ);
      lRIJ  = sqrt(RIJsq);
      lRKJ  = sqrt(RKJsq);
      /* get the cosine relation */
      w = dotprod(RIJ, RKJ)/lRIJ/lRKJ;
      /* test arccos */
      if (fabs(w) > 1.0)
         w = -w;
      theta = acos(w);
      /* calculate dcos/dx */ 
      dtdw = -1/sqrt(1 - w*w);
      dtdRIJ.x = dtdw*(RKJ.x/lRIJ/lRKJ - w*RIJ.x/RIJsq);
      dtdRIJ.y = dtdw*(RKJ.y/lRIJ/lRKJ - w*RIJ.y/RIJsq);
      dtdRIJ.z = dtdw*(RKJ.z/lRIJ/lRKJ - w*RIJ.z/RIJsq);
      dtdRKJ.x = dtdw*(RIJ.x/lRIJ/lRKJ - w*RKJ.x/RKJsq);
      dtdRKJ.y = dtdw*(RIJ.y/lRIJ/lRKJ - w*RKJ.y/RKJsq);
      dtdRKJ.z = dtdw*(RIJ.z/lRIJ/lRKJ - w*RKJ.z/RKJsq);

      /* calculate energy and force factor */
      switch(peftype)
      {  case 1: tt  =  (theta - T[i].t0);
                 V  +=  T[i].Kt * tt * tt;
                 dV  =  T[i].Kt * tt;
                 break;
         case 2: tt  =  (theta - T[i].t0);
                 if (tt < 0.0)
                 {  V  +=  T[i].Kt * tt * tt;
                    dV  =  T[i].Kt * tt;
                 }
                 else 
                    dV = 0.0;
                 break;
        default: dV  =  0.0;
      }

      /* calculate forces */
      F[ai].x -= dtdRIJ.x * dV;
      F[ai].y -= dtdRIJ.y * dV;
      F[ai].z -= dtdRIJ.z * dV;
      F[aj].x += (dtdRIJ.x + dtdRKJ.x) * dV;
      F[aj].y += (dtdRIJ.y + dtdRKJ.y) * dV;
      F[aj].z += (dtdRIJ.z + dtdRKJ.z) * dV;
      F[ak].x -= dtdRKJ.x * dV;
      F[ak].y -= dtdRKJ.y * dV;
      F[ak].z -= dtdRKJ.z * dV;
   } 
   
   tv += cput() - tstart;
   return (V/2);

}  /* End of Vtheta */


/* End of file */
