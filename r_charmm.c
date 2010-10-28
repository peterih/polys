/**********************************************************************
  "r_charmm.c"

  Written by Soren Balling Engelsen, INRA-93, 94.
**********************************************************************/
#include <string.h>
#include "defs.h"
#include "typs.h"
#include "extern.h"
#include "sinput.h"
#include "fileoper.h"
 

/**********************************************************************
        Reads CHARMM format coordinate files.
        =====================================
lines begining with a asterisk '*' are header (comment) lines
1 LINE:  number of atoms (lines following)
* LINE:  atom number (%5d)
         residue number (%5d)
         residue name (x%-4s)
         atom label (x%-4s)
         X coordinate (%10.5lf)
         Y coordinate (%10.5lf)
         Z coordinate (%10.5lf)
         residue name (x%-4s)        ?
         residue number (x%d)        ?
         partial charge (x%12.5lf)
**********************************************************************/
void read_CHARMM(char *fname)

{  FILE    *fp;
   char    line[MAXLEN+1];
   char    resname[4];
   int     i, idum, nat;

   if ((fp = fileopen(fname, "rt")) != NULL)
   {  printf("CHARMM coordinate file %s opened for read-only access\n",fname);

      /* skip comment lines */
      do
      { idum = get_line(fp, line, MAXLEN);
        printf("LINE: %s\n", line);
      } while (line[0] == '*');

      /* read number of atoms */
      idum = sscanf(line, "%d", &nat);

      printf("Number of atoms in residue %s[%d] = %d\n",M.resid,M.nres,nat);

      if (nat > MAXATOM)
      {  printf("Number of atoms (%d) in file '%s' exceeds the ", nat, fname);
         printf("the allowed maximum of (%d)\n", MAXATOM);
         printf("POLYS is quiting this operation\n");
         return;
      }   

      for ( i=0; i<nat; i++)
      {  fscanf(fp, "%d%d%s%s%lf%lf%lf%s%d%lf",
                &idum, &A[i].res, &resname, &A[i].lab, &A[i].pos.x, 
                &A[i].pos.y, &A[i].pos.z, &resname, &idum, &A[i].chg);
         /* CHARMM coordinate files does not contain atomtype   */
         /* information - lets try to estimate atom types       */
         switch (A[i].lab[0])
         {  case 'h': A[i].ztyp =  101; break;
            case 'H': A[i].ztyp =  101; break;
            case 'o': A[i].ztyp =  802; break;
            case 'O': A[i].ztyp =  802; break;
            case 'c': A[i].ztyp =  604; break;
            case 'C': A[i].ztyp =  604; break;
            case 'n': A[i].ztyp =  703; break;
            case 'N': A[i].ztyp =  703; break;
            default : A[i].ztyp =    9; break;
         }

         M.nat++;
      }
      fileclose(fname, fp);

      M.nres = A[M.nat-1].res;

      printf("Number of atoms in %s = %d\n", M.id, M.nat);
      printf("Coordinates for %s succesfully read from file %s\n\n", 
                              M.id,                       fname);
   }
} /* End of read_CHARMM */

/* End of file */
