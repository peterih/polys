/*
 * This file contains the source code for transformations of cartesian 
 * coordinates into internal coordinates.
 * Written by Soren Balling Engelsen, INRA-93.
 */
#include <stdio.h>
#include <math.h>
#include "polys.h"
#include "extern.h"

/*
 * Function angval calculates the angle between the points a-b-c in 
 * cartesian space. Using the formula. 
 *                         rabúrbc
 *    cos(theta) = PI - -------------
 *                       |rab|ú|rbc|
 */
double
angval(Vector3 a, Vector3 b, Vector3 c)
{
   Vector3  vab, vbc;
   double   theta;

   /* get the two bond vectors */
   vab.x = a.x - b.x;
   vab.y = a.y - b.y;
   vab.z = a.z - b.z;
   vbc.x = b.x - c.x;
   vbc.y = b.y - c.y;
   vbc.z = b.z - c.z;
   /* normalize */
   normalize(&vab);
   normalize(&vbc);
   /* get the cosine relation */
   theta = dotprod(vab, vbc);
   /* test arccos */
   if (fabs(theta) > 1.0)
      theta = -theta;
   theta = 180.0 - acos(theta)*RADtoDEG;
   return theta;
}

/*
 * Generates angletabel T[] and counts the number of 1-3 inter-
 * actions in the molecule. This tabel should be used for force fields 
 * which needs to ignore geminal nonbonded interactions.
 */
void
gen_angletab()
{
   register int bi;
   register int bj;
   int          ano;

   printf("Generating angle tabel for %s\n\n", M.id);
 
   ano = 0;                                   /* reset number of angles */
   for ( bi=0; bi<M.nbond-1; bi++)          /* outer loop over bonds */
      for ( bj=bi+1; bj<M.nbond; bj++)      /* inner loop over bonds */
      {  if (B[bi].from == B[bj].from)
         {  T[ano].from = B[bi].to;
            T[ano].over = B[bi].from;
            T[ano++].to = B[bj].to;
         }
         if (B[bi].to   == B[bj].from)
         {  T[ano].from = B[bi].from;
            T[ano].over = B[bi].to;
            T[ano++].to = B[bj].to;
         }
         if (B[bi].from == B[bj].to)
         {  T[ano].from = B[bi].to;
            T[ano].over = B[bi].from;
            T[ano++].to = B[bj].from;
         }
         if (B[bi].to   == B[bj].to)
         {  T[ano].from = B[bi].from;
            T[ano].over = B[bi].to;
            T[ano++].to = B[bj].from;
         }
      }
   M.ntheta = ano;

   printf("  Number of valence angles            : %7d\n\n", M.ntheta);
}

/*
 * Function torval calculates the torsional angle between the points 
 * a-b-c-d in cartesian space. Using the formula. 
 *                     (rab x rbc)ú(rbc x rcd)
 *    cos(phi) = PI - -------------------------
 *                           |rab|ú|rbc|
 */
double
torval(Vector3 a, Vector3 b, Vector3 c, Vector3 d)
{
   Vector3  vab, vbc, vcd, vabc, vbcd, vabcd;
   double   phi, labc, lbcd;

   /* get the three bond vectors */
   vab.x = a.x - b.x;
   vab.y = a.y - b.y;
   vab.z = a.z - b.z;
   vbc.x = b.x - c.x;
   vbc.y = b.y - c.y;
   vbc.z = b.z - c.z;
   vcd.x = c.x - d.x;
   vcd.y = c.y - d.y;
   vcd.z = c.z - d.z;
   /* calculate crossproducts */
   crossprod(vab, vbc, &vabc);
   crossprod(vbc, vcd, &vbcd);
   /* normalize */
   labc = normalize(&vabc);
   lbcd = normalize(&vbcd);
   if (labc*lbcd < RODEPS)
      phi = 0.0;
   else
      /* calculate dotproduct */
      phi = -dotprod(vabc, vbcd);
   /* test arccos */
   if (phi > 1.0)
      phi = 1.0;
   if (phi < -1.0)
      phi = -1.0;
   phi = 180.0 - acos(phi)*RADtoDEG;
   /* determination of sign of phi */
   crossprod(vabc, vbcd, &vabcd);
   normalize(&vabcd);
   if (dotprod(vabcd, vbc) > 0.0)
      phi = -phi;
   return phi;
}

/* Updates the torsional angle tabel P[] */
static void
new_tor(int a, int b, int c, int d, int *tnr) 
{
   P[*tnr].from  = a;
   P[*tnr].over1 = b;
   P[*tnr].over2 = c;
   P[*tnr].to    = d;
   (*tnr)++;
}

/*
 * Generates tabel of torsional angles P[] and counts the number
 * of 1-4 interactions in the molecule. This tabel should be used for force 
 * fields which needs to scale 1-4 nonbonded interactions.
 */
void
gen_tortab()
{
   register int bi;
   register int bj; 
   register int bk; 
   int      tno;

   printf("Generating torsional angle tabel for %s\n\n", M.id);

   tno = 0;                                   /* reset number of torsions */
   for ( bi=0; bi<M.nbond; bi++)            /* loop over central bond */
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
                           new_tor(B[bj].to, B[bi].from,
                                   B[bi].to, B[bk].to, &tno);
                        if (B[bk].to==B[bi].to)
                           new_tor(B[bj].to, B[bi].from,
                                   B[bi].to, B[bk].from, &tno);
                     }
               } 
               /* CASE 2 */
               if (B[bj].to==B[bi].from)
               {  /* connection forwards */
                  for ( bk=0; bk<M.nbond; bk++)
                     if ((bk!=bj) AND (bk!=bi))
                     {  if (B[bk].from==B[bi].to)
                           new_tor(B[bj].from, B[bi].from, 
                                   B[bi].to, B[bk].to, &tno);
                        if (B[bk].to==B[bi].to)
                           new_tor(B[bj].from, B[bi].from,
                                   B[bi].to, B[bk].from, &tno);
                     }
               }
            }

   M.nphi = tno;
   printf("  Number of torsional angles          : %7d\n\n", M.nphi);
}

/* */
void
Add_Element(int *count, int expt, int exa, int index[], int TA[]) 
{
   BOOLEAN new;
   int     n;

   new = TRUE;
   /* check if atom "exa" already has been included in the 
      exclusion list for atom ai                             */
   for (n=index[expt]; n < index[expt+1]; n++)
      if (TA[n] == exa)
         new = FALSE;
   if (new)
   {  TA[index[expt + 1]] = exa;
      index[expt+1] = index[expt+1]+1;
      (*count)++;
   }
}

/*
 * Generates a tabel of the complementary maengde to the nonbonded 
 * interactions.
 */
void
gen_nonbond() 
{
   register int i;
   register int j;
   int          nex = 0;
   int          n14 = 0;
   BOOLEAN      onethree = TRUE;

   printf("Generating thin nonbonded exclusion matrix for %s\n\n", M.id);

   INB[0] = 0;
   I14[0] = 0;
   for (i=0; i < M.nat-1; i++)               /* outer loop over atoms */
   {  INB[i+1] = INB[i];
      I14[i+1] = I14[i];
      for (j=0; j < M.nbond; j++)             /* loop over bonds */
         if ((B[j].from == i) AND (B[j].to > i))
            Add_Element(&nex, i, B[j].to, INB, ANB);
      if (onethree)
         for (j=0; j < M.ntheta; j++)         /* loop over angles */
         {   if ((T[j].from == i) AND (T[j].to > i))
                Add_Element(&nex, i, T[j].to, INB, ANB);
             if ((T[j].to == i)   AND (T[j].from > i))
                Add_Element(&nex, i, T[j].from, INB, ANB);
         }
      for (j=0; j < M.nphi; j++)              /* loop over dihedrals */
      {   if ((P[j].from == i) AND (P[j].to > i))
             Add_Element(&n14, i, P[j].to, I14, A14);
          if ((P[j].to == i)   AND (P[j].from > i))
             Add_Element(&n14, i, P[j].from, I14, A14);
      }
      Qsort(INB[i], INB[i+1]-1, ANB);
      Qsort(I14[i], I14[i+1]-1, A14);
   }

   printf("  Number of nonbonded exclusions      : %7d\n", nex);
   printf("  Number of 1-4 nonbonded interactions: %7d\n\n", n14);
}

/*
 * Creates a list of all internal coordinates on standard output.
 * option 0:    only rotatable bonds
 * option 1:    all except nonbonded exclusions
 * option 2:    all
 */
void
list_intcoor(int option)
{
   register int i;
   int    n, m;
   double l, theta, phi;

   if (option > 0)
   {  printf("Number of bonds: %d\n\n", M.nbond);
      for (i=0; i<M.nbond; i++)
      {  l=dist(A[B[i].from].pos, A[B[i].to].pos);
         printf(" BOND[%5d]: %-8s[%3d] %-8s[%3d] == %5.3lf\n", 
                i, A[B[i].from].lab, A[B[i].from].res, 
                   A[B[i].to].lab,   A[B[i].to].res, l);
      }
      printf("\n\nNumber of valence angles: %d\n\n", M.ntheta);
      for (i=0; i<M.ntheta; i++)
      {  theta=angval(A[T[i].from].pos, A[T[i].over].pos, A[T[i].to].pos); 
         printf(" ANGL[%5d]: %-8s[%3d] %-8s[%3d] %-8s[%3d] == %6.2lf\n", 
                i, A[T[i].from].lab, A[T[i].from].res, A[T[i].over].lab,
                   A[T[i].over].res, A[T[i].to].lab, A[T[i].to].res, theta);
      }
      printf("\n\nNumber of torsional angles: %d\n\n", M.nphi);
      for (i=0; i<M.nphi; i++)
      {  phi=torval(A[P[i].from].pos,  A[P[i].over1].pos,
                    A[P[i].over2].pos, A[P[i].to].pos);
         printf(" TORS[%5d]: %-8s[%3d] %-8s[%3d] %-8s[%3d] %-8s[%3d] == %6.1lf\n", 
                        i, A[P[i].from].lab, A[P[i].from].res, 
                           A[P[i].over1].lab, A[P[i].over1].res,
                           A[P[i].over2].lab, A[P[i].over2].res,
                           A[P[i].to].lab,    A[P[i].to].res, phi);
      }
   }
   printf("\n\nNumber of rotatable bonds: %d\n\n", M.nrot);
   printf("\n\nAtom numbers have been incremented by 1 for Cicada\n");
   for (i=0; i<M.nrot; i++)
   {  phi=torval(A[RT[i].from].pos,  A[RT[i].over1].pos,
                 A[RT[i].over2].pos, A[RT[i].to].pos);
      printf(" ROTB[%2d]:      %-4s - %-4s - %-4s - %-4s     =[%2d]= %6.1lf\n",
                     i, A[RT[i].from].lab, A[RT[i].over1].lab,
                        A[RT[i].over2].lab, A[RT[i].to].lab,
                        abs(RT[i].fold), phi);
      printf("cicada%5i%5i%5i%5i%5i%5i\n",
                        RT[i].from+1, RT[i].over1+1,
                        RT[i].over2+1, RT[i].to+1,
                        1, 0);
   }
   if (option == 2)
   {  printf("\n\nNo nonbonded interaction between:\n\n");
      for (i=0; i < M.nat-1; i++)
         if (INB[i] != INB[i+1])
         {  printf("\n%-9s and ", A[i].lab);
            n = INB[i];
            do
            {  m = 1;
               while ((m<8) AND (n < INB[i+1]))
               {  printf("%-9s", A[ANB[n]].lab);
                  m++;
                  n++;
               }
               printf("\n              ");
            }
            while (n != INB[i+1]);
         }
   }      
   printf("\n\n");
   
   if (option == 2)
   {  printf("\n\n14 nonbonded interaction between:\n\n");
      for (i=0; i < M.nat-1; i++)
         if (I14[i] != I14[i+1])
         {  printf("\n%-9s and ", A[i].lab);
            n = I14[i];
            do
            {  m = 1;
               while ((m<8) AND (n < I14[i+1]))
               {  printf("%-9s", A[A14[n]].lab);
                  m++;
                  n++;
               }
               printf("\n              ");
            }
            while (n != I14[i+1]);
         }
   }      
   printf("\n\n");
}
