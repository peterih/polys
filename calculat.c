/*
 * This file contains the source code for utility calculations
 * Written by Soren Balling Engelsen, INRA-93.
 */
#include <stdio.h>
#include <math.h>
#include "polys.h"
#include "extern.h"

/* Calculates the distance:   i - j */
void
calc_dist(int i, int j)
{
   double rij;

   rij = dist(A[i].pos, A[j].pos);
   printf("Interatomic distance:\n");
   printf("    %s [%d] -- %s [%d]\n", A[i].lab, i, A[j].lab, j);
   printf("          %.3lf\n", rij);
}

/*
 * Calculates the angle:     i - j - k
 * the two distances:        i - j , j - k 
 */
void
calc_angle(int i, int j, int k)
{
   double rij, rjk, aijk;

   rij = dist(A[i].pos, A[j].pos);
   rjk = dist(A[j].pos, A[k].pos);
   aijk = angval(A[i].pos, A[j].pos, A[k].pos); 
   printf("Interatomic angle:\n");
   printf("    %s [%3d] -- %s [%3d] -- %s [%3d]\n", 
               A[i].lab, i, A[j].lab, j, A[k].lab, k);
   printf("             %.3lf        %.3lf\n", rij, rjk);
   printf("                  %.2lf\n", aijk);
}

/*
 * Calculates the torsional angle:     i - j - k - l
 * the two angles:                     i - j - k ,  j - k -l
 * the three distances:                i- j , j - k, k - l 
 */
void
calc_tors(int i, int j, int k, int l)
{
   double rij, rjk, rkl, aijk, ajkl, tijkl;

   rij = dist(A[i].pos, A[j].pos);
   rjk = dist(A[j].pos, A[k].pos);
   rkl = dist(A[k].pos, A[l].pos);
   aijk = angval(A[i].pos, A[j].pos, A[k].pos); 
   ajkl = angval(A[j].pos, A[k].pos, A[l].pos); 
   tijkl = torval(A[i].pos, A[j].pos, A[k].pos, A[l].pos); 
   printf("Interatomic torsional angle:\n");
   printf("    %s [%3d] -- %s [%3d] -- %s [%3d] -- %s [%3d]\n", 
               A[i].lab, i, A[j].lab, j, A[k].lab, k, A[l].lab, l);
   printf("             %.3lf        %.3lf        %.3lf\n", rij, rjk, rkl);
   printf("                  %.2lf       %.2lf\n", aijk, ajkl);
   printf("                         %.2lf\n", tijkl);
}

/* Prints a distance table on standard output. */
void
dist_tabel()
{
   register int i;
   register int j;
   double rij;

   printf("Calculating distance matrix for %s\n\n", M.id);
   for ( i=0; i<M.nat; i++)
       printf("%-8s", A[i].lab);
   printf("\n");
   for ( i=0; i<M.nat; i++)
       printf("========");
   printf("\n");
   for ( i=0; i<M.nat; i++)
   {  for (j=0; j<i; j++)
         printf("        ");
      for ( j=i; j<M.nat; j++)
      {  rij = dist(A[i].pos, A[j].pos);
         printf("%8.3lf", rij);
      }
      printf("\n");
   }
   putchar('\n');
   putchar('\n');
}

double
atomass(int i)
{
   int atnr;

   atnr=(int) (A[i].ztyp / 100);
   switch(atnr) 
   {  case  1: return 1.0079;   /* H  */
      case  2: return 4.0026;   /* He */
      case  3: return 6.941;    /* Li */
      case  4: return 9.0112;   /* Be */
      case  5: return 10.811;   /* B  */
      case  6: return 12.011;   /* C  */
      case  7: return 14.0067;  /* N  */
      case  8: return 15.9994;  /* O  */
      case  9: return 18.9984;  /* F  */
      case 10: return 20.179;   /* Ne */
      case 11: return 22.9898;  /* Na */
      case 12: return 24.305;   /* Mg */
      case 13: return 26.9815;  /* Al */
      case 14: return 28.0855;  /* Si */
      case 15: return 30.9738;  /* P  */
      case 16: return 32.066;   /* S  */
      case 17: return 35.453;   /* Cl */
      case 18: return 39.948;   /* Ar */
      default: printf("BEWARE: Atomic Weight Errornously Reported\n");
   }
}

void
calc_mass()
{
   int    i, atnr;
   double molw;

   molw=0.0;
   for (i=0; i<M.nat; i++)
      molw += atomass(i);

   printf("\nMolecular Weight for %s = %10.3lf Daltons\n", M.id, molw);
}

void
calc_rgyr()
{
   int    i;
   double cmx, cmy, cmz, rgyr, tmass, amass;
   double xx, yy, zz, r2, mr2i;

   cmx=0.0;
   cmy=0.0;
   cmz=0.0;
   tmass=0.0;
   for (i=0; i<M.nat; i++)
   {  amass=atomass(i);
      cmx += A[i].pos.x*amass;
      cmy += A[i].pos.y*amass;
      cmz += A[i].pos.z*amass;
      tmass += amass;
   }
   cmx /= tmass;
   cmy /= tmass;
   cmz /= tmass;
   printf("\nCenter of mass for %s: %15.5lf %15.5lf %15.5lf\n", 
           M.id, cmx, cmy, cmz);
   mr2i=0.0;
   for (i=0; i<M.nat; i++)
   {  amass=atomass(i);
      xx = A[i].pos.x-cmx;
      yy = A[i].pos.y-cmy;
      zz = A[i].pos.z-cmz;
      r2 = xx*xx+yy*yy+zz*zz;
      mr2i += amass*r2;
   }
   rgyr = sqrt(mr2i/tmass);
   printf("\nRadius of gyration for %s = %10.3lf Angstrom\n", M.id, rgyr);

}

void
calc_vol()
{
   int     i;
   double  dx, dy, dz, d2, r2, xi, yi, zi;
   double  gincr, gincr3, boxvol, allvol, molvol;
   Vector3 min, max, gp;

   gincr =0.1;
   gincr3=gincr*gincr*gincr;
   molvol=0.0;

   /* first generate grid */
   max.x = -9999.0;
   max.y = -9999.0;
   max.z = -9999.0;
   min.x =  9999.0;
   min.y =  9999.0;
   min.z =  9999.0;

   for (i=0; i<M.nat; i++)
   {  min.x = MIN(A[i].pos.x,min.x);
      min.y = MIN(A[i].pos.y,min.y);
      min.z = MIN(A[i].pos.z,min.z);
      max.x = MAX(A[i].pos.x,max.x);
      max.y = MAX(A[i].pos.y,max.y);
      max.z = MAX(A[i].pos.z,max.z);
   }
   printf("\nXRANGE   %10.3lf - %10.3lf A\n", min.x, max.x);
   printf("\nYRANGE   %10.3lf - %10.3lf A\n", min.y, max.y);
   printf("\nZRANGE   %10.3lf - %10.3lf A\n", min.z, max.z);
   min.x = min.x - 2.0;
   min.y = min.y - 2.0;
   min.z = min.z - 2.0;
   max.x = max.x + 2.0;
   max.y = max.y + 2.0;
   max.z = max.z + 2.0;
   printf("\nXRANGE   %10.3lf - %10.3lf A\n", min.x, max.x);
   printf("\nYRANGE   %10.3lf - %10.3lf A\n", min.y, max.y);
   printf("\nZRANGE   %10.3lf - %10.3lf A\n", min.z, max.z);
   boxvol = (max.x-min.x)*(max.y-min.y)*(max.z-min.z);
   printf("\nBox Volume = %10.3lf A^3\n", boxvol);

   allvol = 0.0;
   for (xi=min.x; xi<max.x; xi+=gincr)
      for (yi=min.y; yi<max.y; yi+=gincr)
         for (zi=min.z; zi<max.z; zi+=gincr)
         {  allvol += gincr3;
            for (i=0; i<M.nat; i++)
            {  dz = zi - A[i].pos.z;
               if (dz < 1.8) 
               {  dy = yi - A[i].pos.y;
                  dx = xi - A[i].pos.x;
                  d2 = dx*dx + dy*dy + dz*dz;
                  switch(A[i].ztyp/100)
                  {  case 1: r2 = 1.21;   /* r(H) = 1.1 */
                             break;
                     case 6: r2 = 3.24;   /* r(C) = 1.8 */
                             break;
                     case 7: r2 = 4.0;    /* r(N) = 2.0 */
                             break;
                     case 8: r2 = 2.56;   /* r(O) = 1.6 */
                             break;
                     default: printf("No van der Waals radii information\n");
                  }
                  if (d2 < r2)
                  {  molvol+=gincr3;
                     break;
                  }
               }
            }
         }


   allvol = (allvol - boxvol)*100 / boxvol;
   printf("\nMolecular Volume %s = %10.3lf A^3\n", M.id, molvol);
   printf("\nEstimated Integration Error = %10.3lf per cent\n", allvol);
   printf("\nOccupied Volume Fraction = %10.3lf\n", molvol/boxvol);
}
