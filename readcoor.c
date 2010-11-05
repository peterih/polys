/*
 * Code for entering Cartesian coordinates into polys.
 * Written by Soren Balling Engelsen, INRA-93.
 * Last revision: S.B.E. 1997
 */
#include <math.h>
#include <string.h>
#include "polys.h"
#include "extern.h"
 
/* Simple routine which deletes the quotes from the quoted string */
void
DelQuotes(char *str)
{
   register int i;
   char     strs[MAXLEN];

   strcpy(strs,str);
   i = 0;
   do
   {  if (strs[i] != '"')
      {  *str = strs[i]; 
         *str++;
      }
   } while (strs[i++] != '\0');
}

/*
 * Guides the input rutines.
 * AMBER, SPF and CFF formats have not been programmed yet and will
 * default to PF format.
 */
void
ReadCoord(int format, char *fname, BOOLEAN ref)
{
   register int i;

   /* first remove the quotes from the file name */
   DelQuotes(fname);

   printf("Reading Cartesian Coordinates in ");
   switch (format)
   {  case 0     : printf("AMBER format\n");
                   printf("AMBER format NOT implemented yet\n");
                   break;
      case 1     : printf("CHARMM format\n");
                   read_CHARMM(fname); 
                   break;
      case 2     : printf("CFF format\n");
                   printf("CFF format NOT implemented yet\n");
                   break;
      case 3     : printf("CSD format\n");
                   read_CSD(fname, ref); 
                   break;
      case 4     : printf("DAT format\n");
                   read_DAT(fname); 
                   break;
      case 5     : printf("MM3 format\n");
                   read_MM3(fname); 
                   break;
      case 6     : printf("MONOBANK format\n");
                   read_MB(fname, FALSE); 
                   break;
      case 7     : printf("PDB format\n");
                   read_PDB(fname);
                   break;
      case 8     : printf("PF format\n");
                   read_PF(fname, ref); 
                   break;
      case 9     : printf("PLOT format\n");
                   printf("The PLOT format CANNOT be exported!\n");
                   break;
      case 10    : printf("SYBYL format\n");
                   read_SYBYL(fname); 
                   break;
      case 11    : printf("SPF format\n");
                   printf("SPF format NOT implemented yet\n");
                   break;
      case 12    : printf("PIM format\n");
                   printf("PIM format NOT implemented yet\n");
                   break;
      case 13    : printf("XYZ format\n");
                   printf("Tinker XYZ format NOT implemented yet\n");
                   break;
   }
}

/*
 * Reads atomic partial charges from a file.
 * 
 * 1 LINE:  comments
 * 2 LINE:  atom number (%5d)
 * * LINE:  partial charge (%10.5lf)
 */
void ReadCharges(char *fname)
{
   FILE    *fp;
   char    line[MAXLEN+1];
   int     i, idum, nchg, unused;
   double  oldchg, sumchg;

   /* first remove the quotes from the file name */
   DelQuotes(fname);

   if ((fp = fileopen(fname, "rt")) != NULL)
   {  printf("PARTIAL CHARGE file %s opened for read-only access\n",fname);

      /* read comments */
      idum = get_line(fp, line, MAXLEN);

      /* read number of charges */
      unused = fscanf(fp, "%d\n", &nchg);

      if (nchg != M.nat)
      {  printf("WARNING:\n");
       printf("WARNING: Number of charges (%d) does not correspond to\n", nchg);
         printf("WARNING:   number of atoms (%d) in the system\n", M.nat);
         printf("WARNING:\n");
      }
      else
      {  sumchg = 0.0;
         for (i=0; i<nchg; i++)
         {  oldchg = A[i].chg;
            idum = get_line(fp, line, MAXLEN);
            sscanf(line, "%lf", &A[i].chg);
            sumchg += A[i].chg;
         }
      }
       
      fileclose(fname, fp);

      printf("Charges for '%s' succesfully read from file '%s'\n",M.id,fname); 
      neutralize(M.nat);
   }
}
