/**********************************************************************
POLYS MODULE <phimin.c>

  This module contains rutines for minimizing energy in dihedral space.

  Written by Soren Balling Engelsen, INRA-93.
  Last revision: S.B.E. 1997
**********************************************************************/
#include <stdio.h>
#include <math.h>
#include "polys.h"
#include "extern.h"


/**********************************************************************
   Updates the vector of rotatable bonds - rotvec[].
**********************************************************************/
static void new_rotmem(int bnr, int a, int b, int c, int d, int *pnr)

{  static   BOOLEAN status = FALSE;
   static   int nsav = 0;         /* saves the bond number */
   static   double wsav = 0.0;    /* saves the weight priority */
   BOOLEAN  ring = FALSE;
   int      nrot, level;
   double   weight;

   if ((ring == TRUE) AND (nsav == bnr))
      return;
   /* discard bonds which are a part of a ring */
   nrot = FindConn(M.nat, M.nbond, b, c, &level);
   if (nrot == 0)
   {  ring = TRUE;
      return;
   }
   else
      ring = FALSE;

   /* calculate pseudo weight */
   weight=atomass(a)+atomass(b)+atomass(c)+atomass(d); 

   /* if new bond then reinitialize static variables */
   if (nsav != bnr)
   {  wsav = 0.0;
      nsav = bnr;
      if (status == TRUE)
         (*pnr)++;
   }

   /* take the heaviest torsion */
   if (weight >= wsav)
   {  RT[*pnr].from  = a;
      RT[*pnr].over1 = b;
      RT[*pnr].over2 = c;
      RT[*pnr].to    = d;
      RT[*pnr].fold  = nrot;
      wsav = weight;
      status = TRUE;
   }
} /* End of new_rotmem */   


/**********************************************************************
   Generates a vector of rotatable bonds - rotvec[]. Only one
torsional angle will be choosen as a reference and preferably one
without terminal hydrogens. Torsions which are a part of a cyclic
structure will be discarded.
**********************************************************************/
void gen_rotvec()

{  register int bi;
   register int bj; 
   register int bk; 
   int      i, pno;
   double   phi;

   printf("Generating a vector of rotatable bonds\n\n");

   pno = 0;                                /* reset number of torsions */
   for ( bi=0; bi<M.nbond; bi++)            /* loop over central bond */
      /* discard bonds including single bonded atoms */
      if ((Valency(A[B[bi].from]) > 1) AND (Valency(A[B[bi].to]) > 1))
         /* connection backwards */
         for ( bj=0; bj<M.nbond; bj++)
            if (bj!=bi)                       /* not the same bond twice */
            {  /* CASE 1 */
               if (B[bj].from==B[bi].from) 
               {  /* connection forwards */
                  for ( bk=0; bk<M.nbond; bk++)
                     if ((bk!=bj) AND (bk!=bi))
                     {  if (B[bk].from==B[bi].to)
                           new_rotmem(bi, B[bj].to, B[bi].from, 
                                      B[bi].to, B[bk].to, &pno);
                        if (B[bk].to==B[bi].to)
                           new_rotmem(bi, B[bj].to, B[bi].from,
                                      B[bi].to, B[bk].from, &pno);
                     }
               } 
               /* CASE 2 */
               if (B[bj].to==B[bi].from)
               {  /* connection forwards */
                  for ( bk=0; bk<M.nbond; bk++)
                     if ((bk!=bj) AND (bk!=bi))
                     {  if (B[bk].from==B[bi].to)
                           new_rotmem(bi, B[bj].from, B[bi].from,
                                      B[bi].to, B[bk].to, &pno);
                        if (B[bk].to==B[bi].to)
                           new_rotmem(bi, B[bj].from, B[bi].from,
                                      B[bi].to, B[bk].from, &pno);
                     }
               }
            }

   M.nrot = pno+1;
   printf("  Number of rotatable bonds: %5d\n\n", M.nrot);

} /* End of gen_rotvec */


/**********************************************************************
   The potential energy function as a function of phi instead of as
a function of cartesian coordinates.
**********************************************************************/
double func(double x[])

{  register int  i;
   double        Vpot;
   Tors          Pi;

   /* first rotate */
   for (i=0; i<M.nrot; i++)
   {  Pi = RT[i];
      rotbond(Pi, x[i+1]);
   }

   /* calculate only dihedral and nonbonded energy */
   Vpot = potener(2);

   /* printf("%12.7f\n",Vpot); */
   return Vpot;
 
}  /* End of func */

/**********************************************************************
**********************************************************************/
void phimin(BOOLEAN silent, int maxfun, int maxit, double enetol, double grdtol)

{  int    i, iter, j;
   double delta = 0.02;
   double *p, *g, fret, f1, **xi;
   Tors   Pi;

   /* create some vectors */
   p=dvector(1,M.nrot);
   g=dvector(1,M.nrot);
 
   /* generate the phi vector */
   for (i=1; i<=M.nrot; i++)
   {  Pi = RT[i-1];
      p[i] = torval(A[Pi.from].pos, A[Pi.over1].pos,
                    A[Pi.over2].pos, A[Pi.to].pos);
      /* printf("p[%d] == %.2lf\n", i, p[i]); */
   }

   /* calculate the energy and its derivatives */
   fret = potener(0);
   for (i=1; i<=M.nrot; i++)
   {  Pi = RT[i-1];
      rotbond(Pi, p[i]+delta);
      f1 = potener(2);
      rotbond(Pi, p[i]-delta);
      g[i] = (f1 - fret)/delta;
   }

   /* initialize the matrix and go */
   nfun = 0;
   xi=dmatrix(1,M.nrot,1,M.nrot);
   for (i=1;i<=M.nrot;i++)
      for (j=1;j<=M.nrot;j++)
         xi[i][j]=(i == j ? -g[i] : 0.0);
   powell(silent,p,xi,M.nrot,enetol,&iter,&fret,func);

   printf("POWELL: Number of iterations: %10d\n", iter);
   printf("POWELL: Number of function calls: %6d\n", nfun);
   printf("POWELL: Minimum function value: %16.7lf\n\n", fret);
   free_dmatrix(xi,1,M.nrot,1,M.nrot);
   free_dvector(p,1,M.nrot);
   free_dvector(g,1,M.nrot);

} /* End of phimin */

/* End of file */
