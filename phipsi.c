/**********************************************************************
   "phipsi.c"

  This program .....

  Written by Soren Balling Engelsen, INRA-93.
**********************************************************************/
#include <math.h>
#include <string.h>
#include "defs.h"
#include "sinput.h"
#include "fileoper.h"
 
/**********************************************************************
        Reads polys PHIPSI.MAP files
**********************************************************************/
BOOLEAN debug;

#define    MAXP   1500

//// 20030117 pih: removed 'void' declaration of main.
main()

{  register int i;
   FILE *fp;
   int  nch = 0;
   int  no = 0;
   double   phi[MAXP], psi[MAXP], Eabs[MAXP], Emin, phimin, psimin;
   char line[MAXLEN+1], *fname;

   Emin = 1E10;
   fname = getstring("Write name on recorded PHI-PSI map file: ");
   if ((fp = fileopen(fname, "rt")) != NULL)
   {  printf("PHIPSI map file %s opened for read-only access\n\n", fname);

      do 
      {  nch = getline(fp, line, MAXLEN);
         if (line[0]=='-')
            if (line[1]=='-')
            
            {  nch = sscanf(line, "%*s%lf%lf%lf", 
                                      &phi[no], &psi[no], &Eabs[no]);
               printf("[%4d]   PHI = %.2lf  PSI = %.2lf  Eabs = %.4lf\n",
                        no, phi[no], psi[no], Eabs[no]);
               if (Eabs[no] < Emin)
               {  Emin = Eabs[no];
                  phimin = phi[no];
                  psimin = psi[no];
               }
               no++;
            }

      } while (nch !=-1);
   }
   fileclose(fname, fp);

   printf("\nMINIMUM: PHI= %.2lf  PSI= %.2lf  Emin= %.5lf kcal/mol\n\n", 
                           phimin,     psimin,      Emin);
 
   for (i=0; i<no; i++)
   {  if (phi[i] < 0.0)
         phi[i] += 360.0;
      if (psi[i] < 0.0)
         psi[i] += 360.0;
   }

   fname = getstring("Write name on recorded PHI-PSI map file : ");
   if ((fp = fileopen(fname, "wt")) != NULL)
   {  printf("PHIPSI plot UNIRAS file %s opened for write access\n\n", fname);
      fprintf(fp, "MINIMUM: PHI= %.2lf  PSI= %.2lf  Emin= %.5lf kcal/mol\n",
                                 phimin,     psimin,      Emin);
      fprintf(fp, "     phi       psi            Eabs            Erel\n");
     
      for ( i=0; i<no; i++)
      {  if ((Eabs[i] - Emin) > 10.0)
            Eabs[i] = Emin + 10.0;
         fprintf(fp, "%10.1lf%10.1lf%15.2lf%15.2lf\n",
                      phi[i], psi[i], Eabs[i], Eabs[i] - Emin);
      }
   }
   fileclose(fname, fp);

//// 20030117 pih: added 'return' statement in case of correct termination.
   return 0;
} /* End of main */


/* End of file */
