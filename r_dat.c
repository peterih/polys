/**********************************************************************
  "r_dat.c"

  Written by Soren Balling Engelsen, INRA-93, 94.
**********************************************************************/
#include <string.h>
#include "defs.h"
#include "typs.h"
#include "extern.h"
#include "sinput.h"
#include "fileoper.h"
 
/**********************************************************************
    Simple routine which takes a string and replace all occurencies
of comma's with space's.
**********************************************************************/
void CommaSpace(char *str)

{  register int i;

   while (*str++ != '\0')
      if (*str == ',')
         *str = ' ';

}  /* End of CommaSpace */

/**********************************************************************
        Reads DAT format coordinate files.
        ==================================
1 LINE:  title
2 LINE:  orthogonal transformations
3 LINE:  lines per atom
* LINE:  A[i].lab, A[i].pos.x, A[i].pos.y, A[i].pos.z, type, A[i].chg, res
**********************************************************************/
void read_DAT(char *fname)

{  register int i;
   FILE *fp;
   int  nal, nch = 0;
   char line[MAXLEN+1], uplab[8];
   BOOLEAN orthogonal;

   if ((fp = fileopen(fname, "rt")) != NULL)
   {  printf("Old Format data file %s opened for read-only access\n\n", fname);

      /* ask for residue ID: IUPAC */
      M.resid = getstring("Enter IUPAC id of residue: ");

      /* read titel */
      nch = get_line(fp, line, MAXLEN);
      M.titel = strdup(line);

      /* read orthogonal transformations */
      orthogonal = TRUE;
      nch = get_line(fp, line, MAXLEN);
      CommaSpace(line);
      nch = sscanf(line, "%lf%lf%lf%lf%lf%lf", 
             &CELL.a, &CELL.b, &CELL.c, &CELL.alpha, &CELL.beta, &CELL.gamma);
      if ((CELL.a != 1.0) OR (CELL.b != 1.0) OR (CELL.c != 1.0) OR
          (CELL.alpha != 90.0) OR (CELL.beta != 90.0) OR (CELL.gamma != 90.0))
         orthogonal = FALSE;

      /* read number of lines per atom */
      nch = get_line(fp, line, MAXLEN);
      nch = sscanf(line, "%d", &nal);

      /* read the atomic information */
      i = 0;
      do
      {  nch = get_line(fp, line, MAXLEN);
         CommaSpace(line);
         nch = sscanf(line, "%s%lf%lf%lf%d%lf%d", 
               A[i].lab, &A[i].pos.x, &A[i].pos.y, &A[i].pos.z, 
               &A[i].ztyp, &A[i].chg, &A[i].res);
         printf("[%3d][%3d] TYP [%3d] %s\n",
                 i, A[i].res, A[i].ztyp, A[i].lab);  
         switch(A[i].ztyp)
         /*                                        element         D/A  val */
         {  case  1: A[i].ztyp =   604; break;  /* Carbon          N/N   4  */
            case  2: A[i].ztyp =   802; break;  /* Oxygen    HO-   N/N   2  */
            case  3: A[i].ztyp =   101; break;  /* Hydrogen  CH    N/N   1  */
            case -3: A[i].ztyp =   111; break;  /* Hydrogen  HO-   N/N   1  */
            case  4: A[i].ztyp =   812; break;  /* Oxygen    C-O-C Y/N   2  */
            case  5: A[i].ztyp =   703; break;  /* Nitrogen        Y/N   3  */
            case  6: A[i].ztyp =  1602; break;  /* Sulfur          Y/N   2  */
            case  7: A[i].ztyp =   801; break;  /* Oxygen    C=O   Y/N   2  */
            case 18: A[i].ztyp =  1604; break;  /* Sulfur    SO3-  Y/N   4  */
            default: A[i].ztyp =     9; break;  /* dummy           N/N   -  */
         }
         strcpy(uplab, A[i++].lab);
      } while ((strcmp(uplab, "FIN") != 0) AND (strcmp(uplab, "fin") != 0));
      M.nat += (i-1);
      M.nar = (i-1);

      printf("Number of atoms in residue %d = %d\n", M.nres, M.nat);
      fileclose(fname, fp);

      /* if not orthogonal ---> convert */
      if (!orthogonal)
         frac2carte(CELL, FALSE);

      M.nres = A[i].res + 1;

   }

} /* End of read_DAT */

/* End of file */
