/**********************************************************************
POLYS MODULE <r_mono.c>

  Written by Soren Balling Engelsen, INRA-93,94.
  Last revision: S.B.E. 1997
**********************************************************************/
#include <string.h>
//// 20030116  pih: include of stdlib permits use of getenv in read_MB.
#include <stdlib.h>
#include "defs.h"
#include "typs.h"
#include "extern.h"
#include "sinput.h"
#include "fileoper.h"

/**********************************************************************
        Reads MONOBANK coordinate (binary) files.
        ==========================================
1 LINE:  residue id
2 LINE:  litterature reference
3 LINE:  number of atoms (lines following)
* LINE:  atom number (%6d)
         residue number (%3d)
         X coordinate (%10.5lf)
         Y coordinate (%10.5lf)
         Z coordinate (%10.5lf)
         partial charge (%10.5lf)
         atomtype (%d)
         atom label (%s)
**********************************************************************/
void read_NMB(char *fname)

{  FILE    *fp;
   char    line[MAXLEN+1] = "mono/";
   int     i, idum, nat, nres, unused ;

   /* expand filename */
   strcat(line, fname); 
   strcat(line, ".MB"); 
   if ((fp = fileopen(line, "rb")) != NULL)
   {  /* read residue id */
      idum = get_line(fp, line, MAXLEN);
      M.resid = strdup(line);
      /* read titel */
      idum = get_line(fp, line, MAXLEN);
      /* read number of atoms */
      unused = fscanf(fp, "%5d %5d\n", &nat, &nres); /* number of atoms */
      unused = fread((void *) (A+M.nat), sizeof(A[0]), nat, fp);
      fileclose(fname, fp);
      M.nres++;
      M.nat += nat;
      M.nar = nat;
   }
   else 
   {  printf("READ_NMB:\n"); 
      printf("READ_NMB: Could not find monobank member: %s\n", fname);
      printf("READ_NMB:\n"); 
      exit(0);
   }
} /* End of read_NMB */


/**********************************************************************
        Reads MONOBANK coordinate (ascii) files.
        ==========================================
1 LINE:  residue id
2 LINE:  litterature reference
3 LINE:  number of atoms (lines following)
* LINE:  atom number (%6d)
         residue number (%3d)
         X coordinate (%10.5lf)
         Y coordinate (%10.5lf)
         Z coordinate (%10.5lf)
         partial charge (%10.5lf)
         atomtype (%d)
         atom label (%s)
**********************************************************************/
int read_MB(char *fname, BOOLEAN monob)

{  FILE    *fp;
   char    line[MAXLEN+1];
   int     i, idum, nat, nres, unused;
   char    *pardir = getenv("MONOBANK_DIR");

   if (monob==TRUE)    /* expand filename */
   {  strcpy(line,pardir);
      strcat(line,"/");
      strcat(line,fname);
      strcat(line,".x"); 
   }
   else 
   {  DelQuotes(fname);
      strcpy(line,fname);
      printf("MACRO\n"); 
   }

   if ((fp = fileopen(line, "rt")) != NULL)
   {  idum = get_line(fp, line, MAXLEN);
      M.resid = strdup(line);
      /* read titel */
      idum = get_line(fp, line, MAXLEN);
      /* read number of atoms */
      unused = fscanf(fp, "%5d %5d\n", &nat, &nres); /* number of atoms */
      for ( i=M.nat; i<(M.nat+nat); i++)
      {  unused = fscanf(fp, "%5d%d%lf%lf%lf%lf%d%s",
                  &idum, &A[i].res, &A[i].pos.x, &A[i].pos.y, &A[i].pos.z, 
                  &A[i].chg, &A[i].ztyp, A[i].lab);
         A[i].res += (M.nres-1);
      }
      M.nres+=nres;
      fileclose(line, fp);
      M.nat += nat;
      M.nar = nat;
      return(nres);
   }
   else 
   {  printf("READ_MB:\n"); 
      printf("READ_MB: Could not find MONOBANK member: %s\n", line);
      printf("READ_MB:\n"); 
      exit(0);
   }
} /* End of read_MB */

/* End of file */
