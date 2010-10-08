/**********************************************************************
  "rotate.c"

  This module is the management rutine for rotating about bonds

  Written by Soren Balling Engelsen, INRA-95.
**********************************************************************/
#include <stdlib.h>
#include "polys.h"
#include "extern.h"

/**********************************************************************
**********************************************************************/
void Rotate(int from, int over1, int over2, int to, double newphi)

{  register int  i; 
   register int  j;
   double        oldphi;
   Tors          Pi;

   printf("Modifying Torsional Angle: %s[%d] - %s[%d] - %s[%d] - %s[%d]\n\n",
            A[from].lab, from, A[over1].lab, over1, 
            A[over2].lab, over2, A[to].lab, to);
   oldphi=torval(A[from].pos, A[over1].pos, A[over2].pos, A[to].pos);
   printf("Old value: %10.3lf\n", oldphi);
   Pi.from = from;
   Pi.over1= over1;
   Pi.over2= over2;
   Pi.to   = to;
   rotbond(Pi, newphi);
   oldphi=torval(A[from].pos, A[over1].pos, A[over2].pos, A[to].pos);
   printf("New value: %10.3lf\n", oldphi);
   
} /* End of Rotate */
 
/* End of file */
