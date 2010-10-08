/**********************************************************************
  "Minimize.c"

  This module contains the source for a conjugate gradients algorithm.

  Written by Soren Balling Engelsen, INRA-93.
**********************************************************************/
#include "polys.h"
#include "extern.h"

void
Minimize(int space, int maxfun, int maxit, double enetol,
              double grdtol, BOOLEAN silent)
{
   register int  i;

   if (space == 0)   /* cartesian space */ 
   {  if (!silent)
      {  printf("Energy Minimization (Conjugate Gradients) in Cartesian Space\n");
         wline(1, 40);
         printf("maxfun: %10d   maxit : %10d\n", maxfun, maxit);
         printf("enetol: %10.1e   grdtol: %10.1e\n", enetol, grdtol);
         printf("ndim:   %10d\n", 3*M.nat);
         wline(1, 40);
      }
      cartemini(silent, maxit, enetol);
   }
   if (space == 1)   /* dihedral space */
   {  maxit=10;
      enetol=0.1;
      if (!silent)
      {  printf("Energy Minimization (Powell) in Dihedral Space\n\n");
         wline(1, 40);
         printf("maxfun: %10d   maxit : %10d\n", maxfun, maxit);
         printf("enetol: %10.1e   grdtol: %10.1e\n", enetol, grdtol);
         printf("ndim:   %10d\n", M.nrot);
         wline(1, 40);
      }
      phimin(silent, maxfun, maxit, enetol, grdtol);
   }
}  /* End of Minimize */	
/* End of file */
