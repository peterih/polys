/*
 * This file contains the source code for helix duplications
 * Written by Milou Kouwijzer, INRA-95.
 * Last revision: S.B.E. 97
 */
#include <stdio.h>
#include <math.h>
#include "polys.h"
#include "extern.h"

/* Calculate the difference vector */
void
difvec(Vector3 a, Vector3 b, Vector3 *c)
{
   c->x = a.x - b.x;
   c->y = a.y - b.y;
   c->z = a.z - b.z;

}

/*
 * In this routine helices are duplicated to form double or 
 * triple helices.
 * The number of repeats in the single strand as well as the total
 * number of strands should be given
 */
void
doubleh(int nstr, int nreps, double shift)
{
   int      i, j, nmc=0, nsc=0, ncrsc=0, nrps=0;
   int      rpr, sub_from;
   int      g2r, g2nr, offset=0;
   int      M1=0, M2, M3, M4;
   double   na=0.0, rpt=0.0, apr=0.0, apt=0.0, doo=0.0;
   double   rM1M2, rho, Cl, fact, theta, delta;
   double   a, b, c, cosd, sind, costh, cossh, sinsh, zsh;
   Vector3  B, Bpr, Bdpr, C, Cpr, Eksi, Eeta, Edzeta;

   printf("Duplication of Single Stranded Helix\n\n");

   printf("At least 4 repeating units should have been build\n\n");

   for (i=1; i < M.nres; i++)
     {  if (Res[i].status == 0)
	  nmc++;
	if (Res[i].status == 1)
	  nsc++;
	if (Res[i].status == 2)
	  ncrsc++;
      }
   
   for (i=4; i < M.nres; i++)
     {  if (((M.nres - 1) % i) == 0 AND ((nmc % i) == 0) AND 
	    ((nsc % i) == 0) AND ((ncrsc % i) == 0))
	  nrps = i;
      }
   
   if ((nrps == 0) OR ((nreps > 0 ) AND (nreps < 4)))
     { printf ("Number of repeats too low (should be 4 or more)\n");
       printf("No helical calculations will be done!\n\n\n");
       exit(0);
     }
   
   if (nreps == 0) {
     printf("The helix appears to consist of %d repeats\n",nrps);
     printf("\tIf the number of main chain units and the number of side\n");
     printf("\tchain units have a common divisor, the calculated number of\n");
     printf("\trepeats will probably be too high.\n");
     printf("\tAdd 'nreps=<value>' to the commandline to overcome this.\n\n");
     nreps = nrps;
   }
   else {
     if (((M.nres - 1) % nreps) == 0 AND ((nmc % nreps) == 0) AND 
	 ((nsc % nreps) == 0) AND ((ncrsc % nreps) == 0))
       printf("The helix consists of %d repeats\n\n",nreps);
     else {
       printf("The given number of repeats (%d) seems to be incorrect\n",nreps);
       printf("No helical calculations will be done!\n\n\n");
       exit(0);
     }
   }

   printf("First, calculation of helical parameters\n\n");

   rpr = (M.nres-1) / nreps;

/* from Soren's helix.c */

   g2nr = rpr + 1;
   g2r  = rpr + rpr;
   if (Res[g2r].status != 0)
   {  for (i=g2r; i>g2nr; i--)
        if (Res[i].status == 0)
          break;
      offset = g2r - i;
      g2r = i;
   }

   sub_from = Res[g2r].from;

   switch(sub_from)
   {  case 1: M1 = AtomIndex(g2r, "O_1");
              M2 = AtomIndex(g2nr-offset-1, "O_1");
              break;
      case 2: M1 = AtomIndex(g2r, "O_2");
              M2 = AtomIndex(g2nr-offset-1, "O_2");
              break;
     default: printf("\nHelix: REDUCING ERROR\n");
              exit(0);
   }

   doo = dist(A[M1].pos, A[M2].pos);

/* end from Soren's helix.c */

   /* first orient with first atom of second repeating unit on origin,
      first atom of third repeating unit on x-axis and first atom of first
      repeating unit in xy-plane */

   M1 = 0;
   M2 = AtomIndex(rpr+1, A[0].lab);
   M3 = AtomIndex(2*rpr+1, A[0].lab);
   M4 = AtomIndex(3*rpr+1, A[0].lab);

   Orient(1, M.nres - 1, M2, M3, M1);

   difvec(A[M2].pos, A[M1].pos, &Bpr);
   difvec(A[M3].pos, A[M2].pos, &B);
   difvec(A[M4].pos, A[M3].pos, &Bdpr);

   difvec(Bpr, B, &C);
   difvec(B, Bdpr, &Cpr);

   Cl = length(C);

   costh = dotprod(C,Cpr) / (Cl*Cl);
   rho = Cl / (2 - 2*costh);
   theta = acos(costh);
   na = theta*RADtoDEG;
   rpt = 360 / na;
   crossprod(C, Cpr, &Edzeta);
   apr = dotprod(B, Edzeta) / (Cl*Cl*sin(theta));

   apt=rpt*apr/rpr;
   printf("\n---HELICAL-PARAMETERS--------------------------\n");
   printf("Residues per repeat                : %4d\n", rpr);
   printf("Repeats per helical turn        (n): %7.2lf \n", rpt);
   printf("Helical rotation per repeat    (na): %7.2lf \n", na);
   printf("----------------------------------------------\n");
   printf("Helical repeat extension  (O- - -O): %7.2lf A\n", doo);
   printf("Helical repeat advancement      (h): %7.2lf A\n", apr);
   printf("Helical fiber repeat          (n*h): %7.2lf A\n", apt);
   printf("----------------------------------------------\n\n");

   /* Orient with z-axis coincident with the helix axis and the x-axis pointing
      to the first atom of the second repeating unit */

   printf("Next, orientation of the helix to have the helix axis coincident\n");
   printf("with the z-axis\n\n");

   rM1M2 = dist(A[M1].pos, A[M2].pos);

   fact = Cl*Cl*sin(theta);

   uvec(Bpr, B, &Eksi);
   Edzeta.x /= fact;
   Edzeta.y /= fact;
   Edzeta.z /= fact;
   crossprod(Edzeta, Eksi, &Eeta);

   for (i=0; i<M.nat; i++)
     {
       a = dotprod(Eksi, A[i].pos) + rho;
       b = dotprod(Eeta, A[i].pos);
       c = dotprod(Edzeta, A[i].pos);
       A[i].pos.x = a;
       A[i].pos.y = b;
       A[i].pos.z = c;
       if (debug)
       {  printf("%2d Residue: %d Label: %-8s      x: %8.4lf y: %8.4lf z: %8.4lf\n", 
              i, A[i].res, A[i].lab, A[i].pos.x, A[i].pos.y, A[i].pos.z);
       }
     }

   /* Now duplicate the single strand
      If 'shift' is given, the secon strand will be shifted according 
      to that part of the repeat unit */

   if (nstr > 1) {
     printf("Finally, duplication of the single strand to form a %d-fold helix\n\n",nstr);

     if ((shift > 0.01) AND (nstr = 2)) {
       printf("The second string will be shifted %5.2lf times the repeat unit\n\n",shift);
       zsh = shift*apr;
       cossh = cos(-shift*theta);
       sinsh = sin(-shift*theta);
     }
     else {
       if (nstr > 2)
	 printf("shift in combination with nstr > 2 not implemented yet,\n");
       printf("no shift performed\n\n");
     }

     for (j=1; j<nstr; j++)
       {
	 delta = 2*j*M_PI / nstr;
	 cosd = cos(delta);
	 sind = sin(delta);
	 for (i=0; i<M.nat; i++)
	   {  A[M.nat+((j-1)*M.nat)+i] = A[i];
	      A[M.nat+((j-1)*M.nat)+i].pos.x = A[i].pos.x*cosd+A[i].pos.y*sind;
	      A[M.nat+((j-1)*M.nat)+i].pos.y = A[i].pos.y*cosd-A[i].pos.x*sind;
              /* Update residue numbering */
	      A[M.nat+((j-1)*M.nat)+i].res = A[i].res+j*(M.nres-1);
	   }
       }
      

       if (debug)
       { printf("Residue: %d Label: %-8s      x: %8.4lf y: %8.4lf z: %8.4lf\n", 
	    A[M1].res, A[M1].lab, A[M1].pos.x, A[M1].pos.y, A[M1].pos.z); 
         printf("Residue: %d Label: %-8s      x: %8.4lf y: %8.4lf z: %8.4lf\n", 
	    A[M1+M.nat].res, A[M1+M.nat].lab, A[M1+M.nat].pos.x, 
            A[M1+M.nat].pos.y, A[M1+M.nat].pos.z); 
         printf("Residue: %d Label: %-8s      x: %8.4lf y: %8.4lf z: %8.4lf\n", 
	    A[M2].res, A[M2].lab, A[M2].pos.x, A[M2].pos.y, A[M2].pos.z); 
         printf("Residue: %d Label: %-8s      x: %8.4lf y: %8.4lf z: %8.4lf\n", 
	    A[M2+M.nat].res, A[M2+M.nat].lab, A[M2+M.nat].pos.x, 
            A[M2+M.nat].pos.y, A[M2+M.nat].pos.z); 
        printf("Residue: %d Label: %-8s      x: %8.4lf y: %8.4lf z: %8.4lf\n", 
	    A[M3].res, A[M3].lab, A[M3].pos.x, A[M3].pos.y, A[M3].pos.z); 
        printf("Residue: %d Label: %-8s      x: %8.4lf y: %8.4lf z: %8.4lf\n", 
	    A[M3+M.nat].res, A[M3+M.nat].lab, A[M3+M.nat].pos.x, 
            A[M3+M.nat].pos.y, A[M3+M.nat].pos.z); 
        printf("Residue: %d Label: %-8s      x: %8.4lf y: %8.4lf z: %8.4lf\n", 
	    A[M4].res, A[M4].lab, A[M4].pos.x, A[M4].pos.y, A[M4].pos.z); 
        printf("Residue: %d Label: %-8s      x: %8.4lf y: %8.4lf z: %8.4lf\n", 
	    A[M4+M.nat].res, A[M4+M.nat].lab, A[M4+M.nat].pos.x, 
            A[M4+M.nat].pos.y, A[M4+M.nat].pos.z); 
       }
     if ((shift > 0.01) AND (nstr = 2)) 
     {
       for (i=M.nat; i<2*M.nat; i++) 
       {
	 a = A[i].pos.x*cossh + A[i].pos.y*sinsh;
	 b = A[i].pos.y*cossh - A[i].pos.x*sinsh;
	 c = A[i].pos.z + zsh;
	 A[i].pos.x = a;
	 A[i].pos.y = b;
	 A[i].pos.z = c;
       }

       if (debug)
       { printf("\nResidue: %d Label: %-8s      x: %8.4lf y: %8.4lf z: %8.4lf\n", 
	      A[M1].res, A[M1].lab, A[M1].pos.x, A[M1].pos.y, A[M1].pos.z); 
         printf("Residue: %d Label: %-8s      x: %8.4lf y: %8.4lf z: %8.4lf\n", 
	      A[M1+M.nat].res, A[M1+M.nat].lab, A[M1+M.nat].pos.x, 
              A[M1+M.nat].pos.y, A[M1+M.nat].pos.z); 
         printf("Residue: %d Label: %-8s      x: %8.4lf y: %8.4lf z: %8.4lf\n", 
	      A[M2].res, A[M2].lab, A[M2].pos.x, A[M2].pos.y, A[M2].pos.z); 
         printf("Residue: %d Label: %-8s      x: %8.4lf y: %8.4lf z: %8.4lf\n", 
	      A[M2+M.nat].res, A[M2+M.nat].lab, A[M2+M.nat].pos.x, 
              A[M2+M.nat].pos.y, A[M2+M.nat].pos.z); 
         printf("Residue: %d Label: %-8s      x: %8.4lf y: %8.4lf z: %8.4lf\n", 
	      A[M3].res, A[M3].lab, A[M3].pos.x, A[M3].pos.y, A[M3].pos.z); 
         printf("Residue: %d Label: %-8s      x: %8.4lf y: %8.4lf z: %8.4lf\n", 
	      A[M3+M.nat].res, A[M3+M.nat].lab, A[M3+M.nat].pos.x, 
              A[M3+M.nat].pos.y, A[M3+M.nat].pos.z); 
         printf("Residue: %d Label: %-8s      x: %8.4lf y: %8.4lf z: %8.4lf\n", 
	      A[M4].res, A[M4].lab, A[M4].pos.x, A[M4].pos.y, A[M4].pos.z); 
         printf("Residue: %d Label: %-8s      x: %8.4lf y: %8.4lf z: %8.4lf\n", 
	      A[M4+M.nat].res, A[M4+M.nat].lab, A[M4+M.nat].pos.x, 
              A[M4+M.nat].pos.y, A[M4+M.nat].pos.z); 
       }
     }
     
     /* update atom and residue values */
     M.nbond=0;
     M.nat *= nstr;
     M.nres=nstr*(M.nres-1)+1;  
   }
   else
     printf("nstr = %d; no duplication of the single strand will be performed\n\n", nstr);
}
