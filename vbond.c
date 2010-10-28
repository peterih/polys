/**********************************************************************
  "Vbonds.c"

  This file contains the subroutines for calculation of covalent bonded
energy as well as its derivatives vrt. cartesians (F).

  Written by Soren Balling Engelsen, INRA-93.
**********************************************************************/
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "polys.h"
#include "extern.h"

static int peftype;   /* type of bonded potential */

/**********************************************************************
    Function which creates a unique code for each type of covalent bond.
       Ai - Aj
**********************************************************************/
double bondcode(Bond Bi)

{  int  it, jt;

   it = A[Bi.from].ztyp;
   jt = A[Bi.to].ztyp;
   if (it > jt)
     return (it + jt * MAXTINT);
   else
     return (jt + it * MAXTINT); 

} /* End of bondcode */


/**********************************************************************
    Function which creates a unique code for each type of covalent bond.
       Ai - Aj
**********************************************************************/
double paircode(int it, int jt)

{
   if (it > jt)
     return (it + jt * MAXTINT);
   else
     return (jt + it * MAXTINT); 

} /* End of paircode */


/**********************************************************************
**********************************************************************/
void read_bp(int nbond)

{  FILE   *fp;
   int    i, j, nopara, it, jt;
   int    count = 0;
   double bicode, bjcode, Kb, alpha, l0;
   char   line[MAXLEN];
   char   *pardir = getenv("PARAM_DIR");

   strcpy(line,pardir);
   strcat(line,"/BPARAM.PF");
   if ((fp = fileopen(line, "rt")) != NULL)
   {  get_line(fp, line, MAXLEN);
      printf("BONDED POTENTIAL -> %s\n", line);
      get_line(fp, line, MAXLEN);
      sscanf(line, "%d%d", &nopara, &peftype);

      switch(peftype)
      {  case 0: printf("POTENTIAL IGNORED");
                 break;
         case 1: printf("HARMONIC POTENTIAL");
                 break;
         case 2: printf("MORSE POTENTIAL");
                 break;
         case 3: printf("Fourier serie");
         case 4: printf("Inverse power");
        default: printf("POTENTIAL IGNORED");
      }
      putchar('\n');
      putchar('\n');

      for (i=0; i<nopara; i++)
      {  alpha = 0.0;
         get_line(fp, line, MAXLEN);
         switch(peftype)
         {  case 2: sscanf(line, "%d%d%lf%lf%lf", &it, &jt, &Kb, &alpha, &l0);
                    break;
           default: sscanf(line, "%d%d%lf%lf", &it, &jt, &Kb, &l0);
         }
         bicode = paircode(it, jt);
         for (j=0; j< nbond; j++)
         {  bjcode = bondcode(B[j]);
            if (bjcode == bicode)
            { B[j].Kb    = Kb;
              B[j].alpha = alpha;
              B[j].l0    = l0;
              count++;
            }   
         }
      }
      putchar('\n');
      fileclose("BPARAM.PF", fp);
   }
   if (count !=  nbond)
   {  printf("WARNING: MISSING BOND PARAMETERS\n");
      printf("WARNING: Missing parameters to %d bonds\n", nbond - count);
      for (j=0; j< nbond; j++)
         if (B[j].Kb == 0.0)
            printf("%d  %d - %d \n", j , A[B[j].from].ztyp, A[B[j].to].ztyp);
      printf("WARNING: MISSING BOND PARAMETERS\n\n");
   }

} /* End of read_bp */
 

/**********************************************************************
**********************************************************************/
double Vbond(int nbond)

{  register int i;  
   register int j;  
   int      ai, aj;
   double   tstart;
   double   l, ll, V, dV;
   Vector3  RIJ;
 
   tstart = cput();
   V = 0.0;
   if (peftype != 0)
   for (i=0; i<nbond; i++)
   {  ai = B[i].from;
      aj = B[i].to;
      RIJ.x = A[ai].pos.x - A[aj].pos.x;
      RIJ.y = A[ai].pos.y - A[aj].pos.y;
      RIJ.z = A[ai].pos.z - A[aj].pos.z;
      l  = normalize(&RIJ);

      /* calculate the potential */
      switch(peftype)
      {   case 1:  /* the simple harmonic potential */
                   ll = (l - B[i].l0);
                   V += 0.5 * B[i].Kb * ll * ll;
                   dV = B[i].Kb * ll;  
                   break;
          case 2:  /* the Morse potential */
                   ll = (l - B[i].l0);
                   ll = exp(-B[i].alpha*ll);
                   V += B[i].Kb*ll*(ll - 2.0);
                   dV = -2.0*B[i].alpha*B[i].Kb*ll*(ll - 1.0);
                   break;
      }
 
      /* calculate the forces */
      F[ai].x -= RIJ.x * dV;
      F[ai].y -= RIJ.y * dV;
      F[ai].z -= RIJ.z * dV;

      F[aj].x += RIJ.x * dV;
      F[aj].y += RIJ.y * dV;
      F[aj].z += RIJ.z * dV;
   } 
   
   tbo += cput() - tstart;
   return (V);

}  /* End of Vbond */


/* End of file */
