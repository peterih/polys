/* Written by Soren Balling Engelsen, INRA-93,94 */
#include <string.h>
#include "defs.h"
#include "typs.h"
#include "extern.h"
#include "fileoper.h"

/*
 *      Reads POLYS format coordinate files.
 *      ====================================
 * 1 LINE:  residue id
 * 2 LINE:  litterature reference
 * 3 LINE:  orthogonal transformations (A, B, C, alpha, beta, gamma)
 * 4 LINE:  number of atoms (lines following)
 * * LINE:  atom number (%6d)
 *        residue number (%3d)
 *        X coordinate (%10.5lf)
 *        Y coordinate (%10.5lf)
 *        Z coordinate (%10.5lf)
 *        partial charge (%10.5lf)
 *        atomtype (%d)
 *        atom label (%s)
 */
void
read_PF(char *fname, BOOLEAN ref)
{
   FILE    *fp;
   char    line[MAXLEN+1];
   int     i, idum, nat, nbond, nsymop, unused ;
   BOOLEAN orthogonal;

   if ((fp = fileopen(fname, "rt")) != NULL)
   {  printf("POLYS coordinate file '%s' opened for read-only access\n",fname);

      if (ref == FALSE)
      {  /* read residue ID: IUPAC */
         idum = get_line(fp, line, MAXLEN);
         M.resid = strdup(line);

         /* read titel */
         idum = get_line(fp, line, MAXLEN);
         M.titel = strdup(line);
      }
      else
      {  idum = get_line(fp, line, MAXLEN);
         idum = get_line(fp, line, MAXLEN);
         printf("Reading coordinates into the reference set\n");
      }

      /* read orthogonal transformations */
      orthogonal = TRUE;
      idum = fscanf(fp, "%lf%lf%lf%lf%lf%lf\n", 
             &CELL.a, &CELL.b, &CELL.c, &CELL.alpha, &CELL.beta, &CELL.gamma);
      if ((CELL.a != 1.0) OR (CELL.b != 1.0) OR (CELL.c != 1.0) OR
          (CELL.alpha != 90.0) OR (CELL.beta != 90.0) OR (CELL.gamma != 90.0))
         orthogonal = FALSE;

      /* read number of atoms */
      unused = fscanf(fp, "%5d%5d%5d\n", &nat,&nbond,&nsymop);

      if (ref == TRUE)
      {  if (nat != M.nat) 
         {  printf("Number of atoms (%d) in new reference set does not ", nat);
            printf("correspond to number of atoms in main set (%d)\n", M.nat);
            printf("POLYS is quiting this operation\n");
            return;
         }   
         else if (nat > MAXATOMR)
         {  printf("Number of atoms (%d) in file '%s' exceeds the", nat, fname);
            printf("the allowed maximum (%d) in the reference set\n", MAXATOMR);
            printf("POLYS is quiting this operation\n");
            return;
         }   
         else 
            printf("(%d) Atoms are being read into the reference set\n", nat);
      }
      else
         printf("Number of atoms in residue %s[%d] = %d\n",M.resid,M.nres,nat);

      if (nat > MAXATOM)
      {  printf("Number of atoms (%d) in file '%s' exceeds the ", nat, fname);
         printf("the allowed maximum of (%d)\n", MAXATOM);
         printf("POLYS is quiting this operation\n");
         return;
      }   

      /* read the atomic information */
      for ( i=0; i<nat; i++)
      {  idum = get_line(fp, line, MAXLEN);
         if (ref == TRUE)
         {  sscanf(line, "%d%d%lf%lf%lf",
                         &idum, &idum, &Ar[i].x, &Ar[i].y, &Ar[i].z);
         }
         else
         {  sscanf(line, "%5d%d%lf%lf%lf%lf%d%s",
                       &idum,&A[i].res,&A[i].pos.x,&A[i].pos.y,&A[i].pos.z, 
                       &A[i].chg,&A[i].ztyp,A[i].lab);
            M.nat++;
         }
      }

      /* read the bond information */
      for (i=0; i<nbond; i++)
      {  unused = fscanf(fp, "%5d%5d", &B[i].from, &B[i].to);   
         M.nbond++;
      }

      /* read symmetry operations */
/*      for (i=0; i<4; i++)
      {  idum = get_line(fp, line, MAXLEN);
         sscanf(line, "%lf%lf%lf%lf",&TM0[i][0],&TM0[i][1],
                          &TM0[i][2],&TM0[i][3]);
      } 
      for (i=0; i<4; i++)
      {  idum = get_line(fp, line, MAXLEN);
         sscanf(line, "%lf%lf%lf%lf",&TM1[i][0],&TM1[i][1],
                          &TM1[i][2],&TM1[i][3]);
      } 
      TM1[0][3] = TM1[0][3] * CELL.a;
      TM1[1][3] = TM1[1][3] * CELL.b;
      TM1[2][3] = TM1[2][3] * CELL.c;
      for (i=0; i<4; i++)
      {  idum = get_line(fp, line, MAXLEN);
         sscanf(line, "%lf%lf%lf%lf",&TM2[i][0],&TM2[i][1],
                          &TM2[i][2],&TM2[i][3]);
      } 
      TM2[0][3] = TM2[0][3] * CELL.a;
      TM2[1][3] = TM2[1][3] * CELL.b;
      TM2[2][3] = TM2[2][3] * CELL.c;
      for (i=0; i<4; i++)
      {  idum = get_line(fp, line, MAXLEN);
         sscanf(line, "%lf%lf%lf%lf",&TM3[i][0],&TM3[i][1],
                          &TM3[i][2],&TM3[i][3]);
      } 
      TM3[0][3] = TM3[0][3] * CELL.a;
      TM3[1][3] = TM3[1][3] * CELL.b;
      TM3[2][3] = TM3[2][3] * CELL.c;
*/
      fileclose(fname, fp);
   
      if (ref == TRUE)
      {  /* if not orthogonal ---> convert */
         if (!orthogonal)
            frac2carte(CELL, ref);
      }
      else
      {  M.nres = A[i].res + 1;
         /* if not orthogonal ---> convert */
         if (!orthogonal)
            frac2carte(CELL, ref);
         printf("Number of atoms in %s = %d\n", M.id, M.nat);
         printf("Number of bonds in %s = %d\n", M.id, M.nbond);
         printf("Coordinates for %s succesfully read from file %s\n\n", 
                                 M.id,                       fname);
      }
   }
}
