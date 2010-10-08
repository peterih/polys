/**********************************************************************
  "r_pdb.c"

  Written by Soren Balling Engelsen, KVL 95.
**********************************************************************/
#include <string.h>
#include "defs.h"
#include "typs.h"
#include "extern.h"
#include "sinput.h"
#include "fileoper.h"
 
/**********************************************************************
        Reads PDB format coordinate files.
        ===========================================

**********************************************************************/
void read_PDB(char *fname)
{  FILE    *fp;
   char    cdum, line[100], dumtyp[4], cmd[6];
   int     i, nat=0, idum;
   BOOLEAN orthogonal;
   double  fdum;

   if ((fp = fileopen(fname, "rt")) != NULL)
   {  printf("PDB coordinate file %s opened for read-only access\n",fname);

      do
      {  /* get new line from input stream */
         fgets(line, 100, fp);

         if (strncmp("HEADER", line, 6)==0)
         {  printf("%s\n", line);
            M.id = strdup("    ");
            strncpy(M.id, line + 62, 4);
         }

         if (strncmp("COMPND", line, 6)==0)
         {  printf("%s\n", line);
            M.titel = strdup("                        ");
            strncpy(M.titel, line + 10, 25);
         }

         if (strncmp("SOURCE", line, 6)==0)
            printf("%s\n", line);
          

         if (strncmp("AUTHOR", line, 6)==0)
            printf("%s\n", line);
          
         if (strncmp("REVDAT", line, 6)==0)
            printf("%s\n", line);
          
         if (strncmp("JRNL  ", line, 6)==0)
            printf("%s\n", line);

         if (strncmp("HET   ", line, 6)==0)
            printf("%s\n", line);

         if (strncmp("CRYST1", line, 6)==0)
         {  orthogonal = TRUE;
            sscanf(line, "%s %lf %lf %lf %lf %lf %lf\n", cmd,
             &CELL.a, &CELL.b, &CELL.c, &CELL.alpha, &CELL.beta, &CELL.gamma);
            if ((CELL.a != 1.0) OR (CELL.b != 1.0) OR (CELL.c != 1.0) OR
               (CELL.alpha != 90.0) OR (CELL.beta != 90.0) OR 
               (CELL.gamma != 90.0))
               orthogonal = FALSE;
            printf("%s\n", line);
         }

         if (strncmp("ATOM  ", line, 6)==0)
         {  sscanf(line, "%6s%d%*c%c%s%*9c%d%8lf%8lf%8lf%6lf%6lf\n", 
                          &cmd, &idum, &cdum, &A[nat].lab, &A[nat].res,
                          &A[nat].pos.x, &A[nat].pos.y, &A[nat].pos.z,
                          &fdum, &fdum);
            nat++;
         }

         if (strncmp("HETATM", line, 6)==0)
         {  sscanf(line, "%6s%d%*c%c%s%*9c%d%8lf%8lf%8lf%6lf%6lf\n", 
                          &cmd, &idum, &cdum, &A[nat].lab, &A[nat].res,
                          &A[nat].pos.x, &A[nat].pos.y, &A[nat].pos.z,
                          &fdum, &fdum);
            nat++;
         }

         if (nat > MAXATOM)
         {  printf("Number of atoms (%d) in file '%s' exceeds the ",nat,fname);
              printf("the allowed maximum of (%d)\n", MAXATOM);
              printf("POLYS is quiting this operation\n");
              return;
         }

      } while( !feof(fp) );



      for ( i=0; i<nat; i++)
      {  printf("%-10d %-10s %10.5lf %10.5lf %10.5lf\n",
                i, A[i].lab, A[i].pos.x, A[i].pos.y, A[i].pos.z); 

         switch (A[i].lab[0])
         {  case 'H' : if (strcmp(A[i].lab, "HXX")==0) 
                       {  A[i].ztyp =  101; break; }     /* hydrogen   */

            case 'C' : if (strcmp(A[i].lab, "C")==0)
                       {  A[i].ztyp =  604; break; }     /* sp3 carbon */
                       if (strcmp(A[i].lab, "CA")==0)
                       {  A[i].ztyp =  604; break; }     /* sp3 carbon */
                       if (strcmp(A[i].lab, "CB")==0)
                       {  A[i].ztyp =  604; break; }     /* sp3 carbon */
                       if (strcmp(A[i].lab, "CD")==0)
                       {  A[i].ztyp =  604; break; }     /* sp3 carbon */
                       if (strcmp(A[i].lab, "CD1")==0)
                       {  A[i].ztyp =  604; break; }     /* sp3 carbon */
                       if (strcmp(A[i].lab, "CD2")==0)
                       {  A[i].ztyp =  604; break; }     /* sp3 carbon */
                       if (strcmp(A[i].lab, "CE")==0)
                       {  A[i].ztyp =  604; break; }     /* sp3 carbon */
                       if (strcmp(A[i].lab, "CE1")==0)
                       {  A[i].ztyp =  604; break; }     /* sp3 carbon */
                       if (strcmp(A[i].lab, "CE2")==0)
                       {  A[i].ztyp =  604; break; }     /* sp3 carbon */
                       if (strcmp(A[i].lab, "CE3")==0)
                       {  A[i].ztyp =  604; break; }     /* sp3 carbon */
                       if (strcmp(A[i].lab, "CG")==0)
                       {  A[i].ztyp =  604; break; }     /* sp3 carbon */
                       if (strcmp(A[i].lab, "CG1")==0)
                       {  A[i].ztyp =  604; break; }     /* sp3 carbon */
                       if (strcmp(A[i].lab, "CG2")==0)
                       {  A[i].ztyp =  604; break; }     /* sp3 carbon */
                       if (strcmp(A[i].lab, "CH2")==0)
                       {  A[i].ztyp =  604; break; }     /* sp3 carbon */
                       if (strcmp(A[i].lab, "CZ")==0)
                       {  A[i].ztyp =  604; break; }     /* sp3 carbon */
                       if (strcmp(A[i].lab, "CZ2")==0)
                       {  A[i].ztyp =  604; break; }     /* sp3 carbon */
                       if (strcmp(A[i].lab, "CZ3")==0)
                       {  A[i].ztyp =  604; break; }     /* sp3 carbon */
                       else 
                       {  A[i].ztyp =  604; 
                          printf("Atom #%d with label %s assigned as standard sp3 carbon\n", i, A[i].lab);
                          break; }     /* sp3 oxygen */


            case 'O' : if (strcmp(A[i].lab, "O")==0)
                       {  A[i].ztyp =  802; break; }     /* sp3 oxygen */
                       if (strcmp(A[i].lab, "OE1")==0)
                       {  A[i].ztyp =  802; break; }     /* sp3 oxygen */
                       if (strcmp(A[i].lab, "OE2")==0)
                       {  A[i].ztyp =  802; break; }     /* sp3 oxygen */
                       if (strcmp(A[i].lab, "OD1")==0)
                       {  A[i].ztyp =  802; break; }     /* sp3 oxygen */
                       if (strcmp(A[i].lab, "OD2")==0)
                       {  A[i].ztyp =  802; break; }     /* sp3 oxygen */
                       if (strcmp(A[i].lab, "OG")==0)
                       {  A[i].ztyp =  802; break; }     /* sp3 oxygen */
                       if (strcmp(A[i].lab, "OG1")==0)
                       {  A[i].ztyp =  802; break; }     /* sp3 oxygen */
                       if (strcmp(A[i].lab, "OH")==0)
                       {  A[i].ztyp =  802; break; }     /* sp3 oxygen */
                       if (strcmp(A[i].lab, "OXT")==0)
                       {  A[i].ztyp =  802; break; }     /* sp3 oxygen */
                       else 
                       {  A[i].ztyp =  802; 
                          printf("Atom #%d with label %s assigned as standard sp3 oxygen\n", i, A[i].lab);
                          break; }     /* sp3 oxygen */

            case 'N' : if (strcmp(A[i].lab, "N")==0)
                       {  A[i].ztyp =  703; break; }     /* sp3 nitrogen */
                       if (strcmp(A[i].lab, "ND1")==0)
                       {  A[i].ztyp =  703; break; }     /* sp3 nitrogen */
                       if (strcmp(A[i].lab, "ND2")==0)
                       {  A[i].ztyp =  703; break; }     /* sp3 nitrogen */
                       if (strcmp(A[i].lab, "NE")==0)
                       {  A[i].ztyp =  703; break; }     /* sp3 nitrogen */
                       if (strcmp(A[i].lab, "NE1")==0)
                       {  A[i].ztyp =  703; break; }     /* sp3 nitrogen */
                       if (strcmp(A[i].lab, "NE2")==0)
                       {  A[i].ztyp =  703; break; }     /* sp3 nitrogen */
                       if (strcmp(A[i].lab, "NH1")==0)
                       {  A[i].ztyp =  703; break; }     /* sp3 nitrogen */
                       if (strcmp(A[i].lab, "NH2")==0)
                       {  A[i].ztyp =  703; break; }     /* sp3 nitrogen */
                       if (strcmp(A[i].lab, "NZ")==0)
                       {  A[i].ztyp =  703; break; }     /* sp3 nitrogen */
                       else 
                       {  A[i].ztyp =  703; 
                          printf("Atom #%d with label %s assigned as standard sp3 nitrogen\n", i, A[i].lab);
                          break; }     /* sp3 oxygen */

            case 'P' : if (strcmp(A[i].lab, "PXX")==0)
                       {  A[i].ztyp = 1504; break; }     /* default phosphor */

            case 'S' : if (strcmp(A[i].lab, "SD")==0)
                       {  A[i].ztyp =  1604; break; }    /* sulphate */

            case 'M' : if (strcmp(A[i].lab, "MN")==0)
                       {  A[i].ztyp =  2500; break; }    /* sulphate */

            default :  A[i].ztyp =    9;
                       printf("WARNING: Atomtype %s unknown!\n",A[i].lab);
         }
         M.nat++;
      }

      fileclose(fname, fp);

      M.nres = A[i].res + 1;

      printf("Number of atoms in %s = %d\n", M.id, M.nat);
      printf("Coordinates for %s succesfully read from file %s\n\n", 
                              M.id,                       fname);
  }
} /* End of read_PDB */

/* End of file */
