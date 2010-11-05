/*
 * This module contains the bookeeping of the potential energy 
 * calculations as well as the "tester" subroutine for checking the 
 * implementation of analytical derivatives needed for the minimization
 * algorithms.
 * Written by Soren Balling Engelsen, INRA-93.
 */
#include <stdio.h>
#include <math.h>
#include "polys.h"
#include "extern.h"

/*
 * Calculates the root-mean-square product of a vector with 
 * "no" elements of the type Vector3
 */
double
RMSE(int no, Vector3 V[])
{
   int     i;
   double  rmse = 0.0;

   for (i=0; i<no; i++)
      rmse += dotprod(V[i], V[i]);

   rmse /= no*3;
   return (sqrt(rmse));
}

/* Calculates the potential energy */
double
potener(int option)
{ 
   int    i;
   double Vs, Vvdw, Vc = 0.0;

   /* initialize variables */
   for (i=0; i < M.nat; i++)
   {   F[i].x = 0.0;
       F[i].y = 0.0;
       F[i].z = 0.0;
   } 
   M.Epot = 0.0;

   /* calculate energy distributions */
   if (option == 0)   M.Eb   = Vbond(M.nbond);
   if (option <= 1)   M.Et   = Vtheta(M.ntheta);
   if (option <= 2)   M.Ep   = Vphi(M.nphi);
   if (option <= 3)   M.Econ = Vcon(M.nphi, M.ncon);
   if (option <= 4)
   {                  M.Enbs = Vnonbond(1, M.nat, &Vvdw, &Vc);
                      M.Evdw = Vvdw;
                      M.Ec   = Vc;
   }

   M.Epot = M.Eb + M.Et + M.Ep + M.Econ + M.Evdw + M.Ec;

   nfun++;
   return(M.Epot);
}

void
tester()
{
   int      i;
   double   V0, V1, dVnum, dVana, dVdif;
   double   delta = RODEPS;
   double   sqsum = 0.0;
   Vector3  F0[MAXATOM];

   printf("Analytical and numerical gradients and their difference:\n\n");
   printf("%20s%20s%20s\n", "Numerical", "Analytical", "Difference");
   printf("%60s%15s\n", "Absolute", "Relative");
   
   V0 = potener(0);
   for (i=0; i<M.nat; i++)
   {  printf("%s[%d]\n", A[i].lab, A[i].res);
      A[i].pos.x += delta;
      V1 = potener(0);
      A[i].pos.x -= delta;
      dVana = -F[i].x;
      dVnum = (V1 - V0)/delta;
      dVdif = dVnum - dVana;
      printf("X:  %20.10f%20.10f%15.2e%15.2e\n",dVnum,dVana,dVdif,dVdif/dVnum);
      sqsum += dVdif*dVdif;
      A[i].pos.y += delta;
      V1 = potener(0);
      A[i].pos.y -= delta;
      dVana = -F[i].y;
      dVnum = (V1 - V0)/delta;
      dVdif = dVnum - dVana;
      printf("Y:  %20.10f%20.10f%15.2e%15.2e\n",dVnum,dVana,dVdif,dVdif/dVnum);
      sqsum += dVdif*dVdif;
      A[i].pos.z += delta;
      V1 = potener(0);
      A[i].pos.z -= delta;
      dVana = -F[i].z;
      dVnum = (V1 - V0)/delta;
      dVdif = dVnum - dVana;
      printf("Z:  %20.10f%20.10f%15.2e%15.2e\n",dVnum,dVana,dVdif,dVdif/dVnum);
      sqsum += dVdif*dVdif;
   }
   sqsum /= M.nat*3;
   printf("\n\nThe RootMeanSquared residue sum: %.2e\n\n", sqrt(sqsum));
}

/*
 * Calculates the potential energy and prints the energy 
 * distributions on standard output
 */
void
printener(BOOLEAN deriv)
{ 
   potener(0);
   printf("Potential Energy Contributions for %s: \n", M.id);
   putchar('\n');
   printf("   %-40s %15s\n", "TYPE", "KCAL/MOL");
   putchar('\n');
   printf("   %-40s %15.5f\n", "Covalent bonded energy", M.Eb);
   printf("   %-40s %15.5f\n", "Valence angle energy", M.Et);
   printf("   %-40s %15.5f\n", "Torsional energy", M.Ep);
   printf("   %-40s %15.5f\n", "Out-of-plane energy", M.Eo);
   printf("   %-40s %15.5f\n", "Van der Waals energy", M.Evdw);
   printf("   %-40s %15.5f\n", "Coulombic energy", M.Ec);
   printf("   %-40s %15.5f\n", "Shifted nonbonded energy", M.Enbs);
   printf("   %-40s %15.5f\n", "Constrain energy", M.Econ);
   putchar('\n');
   printf("   %-40s %15.5f\n", "Total potential energy", M.Epot);
   putchar('\n');
   printf("   %-40s %15.8f\n", "RMS-gradient", RMSE(3*M.nat, F));
   putchar('\n');
   if (deriv)
      tester();
}
