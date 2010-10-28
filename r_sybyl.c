/**********************************************************************
  "r_sybyl.c"

  Written by Soren Balling Engelsen, INRA-93, 94.
**********************************************************************/
#include <string.h>
#include "defs.h"
#include "typs.h"
#include "extern.h"
#include "sinput.h"
#include "fileoper.h"
 

/**********************************************************************
        Reads SYBYL "mol2" format coordinate files.
        ===========================================
find @<TRIPOS>MOLECULE
1* LINE:  (%45s)                         titel
2* LINE:  (%5d %5d %5d %5d %5d)          M.nat, M.nbond, 1, 0, 1
find @<TRIPOS>ATOM
3* LINE:  (%7d %-8s%10.4lf%10.4lf%10.4lf %-4s %6d <1>  %12.4lf)
          atom number
          atom label
          X coordinate (%10.4lf)
          Y coordinate (%10.4lf)
          Z coordinate (%10.4lf)
          atom type (%-4s)
          residue number (%6d)
          fractional chrage (%12.4lf)
**********************************************************************/
void read_SYBYL(char *fname)

{  FILE    *fp;
   char    line[MAXLEN+1], styp[4], dumtyp[30];
   int     i, nat, idum;

   if ((fp = fileopen(fname, "rt")) != NULL)
   {  printf("SYBYL coordinate file %s opened for read-only access\n",fname);

      /* find @<TRIPOS>MOLECULE */
      while (line[0] != '@')
        idum = get_line(fp, line, MAXLEN);

      /* get titel and number of atoms */
      idum = get_line(fp, line, MAXLEN);
      M.titel = strdup(line);
      idum = get_line(fp, line, MAXLEN);
      idum = sscanf(line,"%5d", &nat);
 
      /* find @<TRIPOS>ATOM */
      while (line[0] != '@')
        idum = get_line(fp, line, MAXLEN);

      printf("Number of atoms in residue %s[%d] = %d\n",M.resid,M.nres,nat);

      if (nat > MAXATOM)
      {  printf("Number of atoms (%d) in file '%s' exceeds the ", nat, fname);
         printf("the allowed maximum of (%d)\n", MAXATOM);
         printf("POLYS is quiting this operation\n");
         return;
      }   

      for ( i=0; i<nat; i++)
      {  idum = get_line(fp, line, MAXLEN);
         sscanf(line, "%d %s %lf %lf %lf %s %d %s %lf",
           &idum, &A[i].lab, &A[i].pos.x, &A[i].pos.y, &A[i].pos.z, 
           &styp, &A[i].res, &dumtyp, &A[i].chg);
         printf("%-10d %10s %10.5lf %10.5lf %10.5lf\n",
                i, A[i].lab, A[i].pos.x, A[i].pos.y, A[i].pos.z);

         switch (styp[0])
         {  case 'H' : A[i].ztyp =  101; break;       /* hydrogen   */
            case 'C' : if (styp[2] == '3')
                       {  A[i].ztyp =  604; break;    /* sp3 carbon */
                       }
                       else if (styp[2] == 'C')
                       {  A[i].ztyp =  604; break;    /* PIM sp3 carbon */
                       }
                       else if (styp[2] == 'A')
                       {  A[i].ztyp =  614; break;    /* PIM a-anomeric C */
                       }
                       else if (styp[2] == 'B')
                       {  A[i].ztyp =  614; break;    /* PIM b-anomeric C */
                       }
                       else if (styp[2] == '2')
                       {  A[i].ztyp =  603; break;    /* sp2 carbon */
                       }
                       else if (styp[2] == '1')
                       {  A[i].ztyp =  602; break;    /* sp  carbon */
                       }
                       else if (styp[1] == 'l')
                       {  A[i].ztyp = 1701; break;    /* sp  carbon */
                       }
                       else 
                          A[i].ztyp =  604;           /* default carbon */
            case 'N' : if (styp[2] == '3')
                       {  A[i].ztyp =  703; break;    /* sp3 nitrogen */
                       }
                       else if (styp[2] == 'a')
                       {  A[i].ztyp =  703; break;    /* amide nitrogen */
                       }
                       else if (styp[2] == '2')
                       {  A[i].ztyp =  702; break;    /* sp2 nitrogen */
                       }
                       else if (styp[2] == '1')
                       {  A[i].ztyp =  701; break;    /* sp  nitrogen */
                       }
                       else 
                          A[i].ztyp =  703;           /* default nitrogen */
            case 'O' : if (styp[2] == '3')
                       {  A[i].ztyp =  802; break;    /* sp3 oxygen */
                       }
                       else if (styp[2] == 'A')
                       {  A[i].ztyp =  812; break;    /* PIM a-anomeric O */
                       }
                       else if (styp[2] == 'B')
                       {  A[i].ztyp =  812; break;    /* PIM b-anomeric O */
                       }
                       else if (styp[2] == 'R')
                       {  A[i].ztyp =  812; break;    /* PIM sugar O-ring */
                       }
                       else if (styp[2] == '2')
                       {  A[i].ztyp =  801; break;    /* sp2 oxygen */
                       }
                       else 
                          A[i].ztyp =  802;           /* default oxygen */
            case 'P' : 
                          A[i].ztyp = 1504;           /* default phosphor */
            case 'S' : if (styp[2] == 'o')
                       {  A[i].ztyp =  1604; break;   /* sulphate */
                       }
                       else if (styp[2] == '3')
                       {  A[i].ztyp =  1604; break;    /* -SO3 sulphur */
                       }
                       else
                          A[i].ztyp = 1602;           /* default sulphur */

            default :     A[i].ztyp =    9;
                      printf("WARNING: Atomtype %s unknown!\n",styp);
         }
         M.nat++;
      }
      fileclose(fname, fp);

      M.nres = A[i].res + 1;

      printf("Number of atoms in %s = %d\n", M.id, M.nat);
      printf("Coordinates for %s succesfully read from file %s\n\n", 
                              M.id,                       fname);
  }
} /* End of read_SYBYL */

/* End of file */
