/*
 * Written by Soren Balling Engelsen, INRA-93.
 */
#include <math.h>
#include <string.h>
#include "defs.h"
#include "typs.h"
#include "sinput.h"
#include "fileoper.h"
 
BOOLEAN debug;
/**********************************************************************
        Translates cartesian coordinates.
**********************************************************************/
#define    MAXA   1500

/* 20030117 pih: removed 'void' declaration of main */
int
main()
{
   register int i;
   FILE   *fp;
   char   line1[80];
   char   line2[80];
   char   line3[80];
   char   *fname;
   int    nch, idum, nat, nres, ntrans;
   double xmin, xmax, xav, xtr; 
   double ymin, ymax, yav, ytr; 
   double zmin, zmax, zav, ztr; 
   atom   A[MAXA];

   xmin = 1E10; ymin = 1E10; zmin = 1E10;
   xmax =-1E10; ymax =-1E10; zmax =-1E10;
   xav = 0.0;
   yav = 0.0;
   zav = 0.0;
   fname = getstring("Write name on polys coordinate file: ");
   if ((fp = fileopen(fname, "rt")) != NULL)
   {  printf("POLYS coordinate file %s opened for read-only access\n\n", fname);
      nch = get_line(fp, line1, MAXLEN);  /* ID */
      nch = get_line(fp, line2, MAXLEN);  /* title */
      nch = get_line(fp, line3, MAXLEN);  /* crystal */
      fscanf(fp, "%d", &nat);
      printf("NAT: %d\n", nat);
      for (i=0; i<nat; i++)
      {  fscanf(fp, "%5d%d%lf%lf%lf%lf%d%s",
                &idum, &A[i].res, &A[i].pos.x, &A[i].pos.y, &A[i].pos.z,
                &A[i].chg, &A[i].ztyp, &A[i].lab);
         if (A[i].pos.x < xmin)
            xmin = A[i].pos.x;
         if (A[i].pos.y < ymin)
            ymin = A[i].pos.y;
         if (A[i].pos.z < zmin)
            zmin = A[i].pos.z;
         if (A[i].pos.x > xmax)
            xmax = A[i].pos.x;
         if (A[i].pos.y > ymax)
            ymax = A[i].pos.y;
         if (A[i].pos.z > zmax)
            zmax = A[i].pos.z;
         xav += A[i].pos.x;
         yav += A[i].pos.y;
         zav += A[i].pos.z;
      }
      xav /= nat;
      yav /= nat;
      zav /= nat;
   }
   fileclose(fname, fp);

   printf("\nCOORDINATE        X           Y            Z\n");
   printf("minimum %12.5lf %12.5lf %12.5lf\n", xmin, ymin, zmin);
   printf("maximum %12.5lf %12.5lf %12.5lf\n", xmax, ymax, zmax);
   printf("average %12.5lf %12.5lf %12.5lf\n", xav,  yav,  zav);

   ntrans=0;
   printf("Write desired coordinate translations (x y z<cr>):\n");
   printf("%d:    0.0     0.0    0.0\n", ntrans);
   do
   {  ntrans++;
      printf("%d: ", ntrans);
      nch = scanf("%lf%lf%lf", &xtr, &ytr, &ztr);
      printf("%d:  %10.2lf%10.2lf%10.2lf\n", ntrans, xtr, ytr, ztr);
      if (nch == 3)
      {  for (i=0; i<nat; i++)
         {  A[i+ntrans*nat] = A[i];
            A[i+ntrans*nat].pos.x = A[i].pos.x + xtr;
            A[i+ntrans*nat].pos.y = A[i].pos.y + ytr;
            A[i+ntrans*nat].pos.z = A[i].pos.z + ztr;
         }
      }
   } while (nch == 3);
   nat = (ntrans)*nat;

   fname = getstring("Write name on recorded new coordinate file: ");
   if ((fp = fileopen(fname, "wt")) != NULL)
   {  printf("POLYS coordinate file %s opened for write access\n\n", fname);
      fprintf(fp, "%s\n", line1);
      fprintf(fp, "%s\n", line2);
      fprintf(fp, "%s\n", line3);
      fprintf(fp, "%5d\n", nat);
      for ( i=0; i<nat; i++)
         fprintf(fp, "%5d %3d %15.5f %15.5f %15.5f %7.4f %6d %s\n",
                     i, A[i].res, A[i].pos.x, A[i].pos.y, A[i].pos.z,
                        A[i].chg, A[i].ztyp, A[i].lab);
   }
   fileclose(fname, fp);
// 20030117 pih: added a 'return' statement in case of correct termination.
   return 0; 
} /* End of main */
