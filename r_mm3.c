/* Written by Soren Balling Engelsen, INRA-93,94 */
#include <string.h>
#include "defs.h"
#include "typs.h"
#include "extern.h"
#include "sinput.h"
#include "fileoper.h"

/*
 *      Reads MM3 format coordinate files.
 *      =====================================
 * 1* LINE:  (6A10,I1,I4,I2,I1,2I2,I3,F5.0) Name card
 * 1a* LINE: (60L1,2I1,7I2,2I1,I2)          Conjugated phi-system
 * 2* LINE:  (I1,I4,5X,F10.5,9I5,5X,2I5)    Lists and options
 * 5* LINE:  (3F10.5, I5)                   Coordinates
 *       X coordinate (%10.5lf)
 *       Y coordinate (%10.5lf)
 *       Z coordinate (%10.5lf)
 *       atom type letter (%c)           Imberty
 *       atom type (%2d)
 *       atom label (%s)                 Imberty
 */
void
read_MM3(char *fname)
{
   FILE    *fp;
   char    line[MAXLEN+1], titel1[MAXLEN+1], titel2[MAXLEN+1];
   char    cha;
   int     i, idum, method, nat, kfixtyp, ncon, nattch, nls, attyp;
   int    unused ;
   double  attch;

   if ((fp = fileopen(fname, "rt")) != NULL)
   {  printf("MM3 coordinate file %s opened for read-only access\n",fname);

      /* get title and number of atoms from line 1* */
      idum = get_line(fp, line, MAXLEN);
      idum = sscanf(line,"%60c%1d%4d", titel1, &method, &nat);
      idum = sscanf(titel1, "%s", titel2);
      M.titel = strdup(titel2);

      /* if pi system */
      if (method==1)
         idum = get_line(fp, line, MAXLEN);

      /* get number of lines with connctivity from line 2* */
      idum = get_line(fp, line, MAXLEN);
      idum = sscanf(line, "%1d%4d%22c%5d", &kfixtyp,&ncon,titel1,&nattch);
      nls = ncon + (int) nattch/8;
/* Old criterion did not work properly!
      if ((fmod(nattch, 8.0)) > 0.0) 
	nls++; */
      if (nattch - 8*((int) nattch/8) > 0)
	nls++;

      printf("Number of atoms in residue %s[%d] = %d\n",M.resid,M.nres,nat);
      printf("Number of lines to skip = %d + %d/8 = %d\n", ncon, nattch, nls); 

      for ( i=0; i<nls; i++)
         idum = get_line(fp, line, MAXLEN);

      for ( i=0; i<nat; i++)
      {  unused = fscanf(fp, "%lf%lf%lf %c %2d%s",
                  &A[i].pos.x,&A[i].pos.y,&A[i].pos.z,&cha,&attyp,A[i].lab);

         switch (attyp)
         {  case  1: A[i].ztyp =  604; break;    /* sp3 */
            case  2: A[i].ztyp =  603; break;    /* sp2 - alkene */
            case  3: A[i].ztyp =  603; break;    /* sp2 - carbonyl */
            case  4: A[i].ztyp =  602; break;    /* sp  - alkyne */
            case  5: A[i].ztyp =  101; break;    /* except on N, O, S */
            case  6: A[i].ztyp =  802; break;    /* C-O-H C-O-C O-O;
                                                    beware: both 802 and 812 
                                                    are possible here! */
            case  7: A[i].ztyp =  801; break;    /* C=O carbonyl */
            case  8: A[i].ztyp =  713; break;    /* sp3 */
            case  9: A[i].ztyp =  703; break;    /* sp2 */
            case 10: A[i].ztyp =  701; break;    /* sp  */
            case 15: A[i].ztyp = 1602; break;    /* -S- sulfide */
            case 17: A[i].ztyp = 1603; break;    /* sulfur oxide */
            case 21: A[i].ztyp =  111; break;    /* -OH alcohol */
            case 23: A[i].ztyp =  131; break;    /* NH amine */
            case 24: A[i].ztyp =  111; break;    /* COOH  carboxyl */
            case 28: A[i].ztyp =  121; break;    /* NH amide */
            case 42: A[i].ztyp = 1602; break;    /* sp2 thiophene */
            case 75: A[i].ztyp =  822; break;    /* O-H, O-C carboxyl */
            case 77: A[i].ztyp =  801; break;    /* O=C- carboxyl acid */
            case 78: A[i].ztyp =  801; break;    /* O=C- carboxyl ester */
            case 79: A[i].ztyp =  801; break;    /* O=C-N< amide */
            default : A[i].ztyp =   9;
                      printf("WARNING: Atomtype not specified properly!\n");
         }

         M.nat++;
      }
      fileclose(fname, fp);

      M.nres = A[i].res + 1;

      printf("Number of atoms in %s = %d\n", M.id, M.nat);
      printf("Coordinates for %s succesfully read from file %s\n\n", 
                              M.id,                       fname);
   }
}
