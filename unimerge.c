/**********************************************************************
  "mergeuni.c"

  This program .....

  Written by Soren Balling Engelsen, INRA-93.
**********************************************************************/
#include <math.h>
#include <string.h>
#include "defs.h"
#include "sinput.h"
#include "fileoper.h"
 
/**********************************************************************
        Reads uniras *.uni files
**********************************************************************/
BOOLEAN  debug;

#define    MAXP   1500

//// 20030117 pih: removed 'void' declaration of main.
int main()

{  register int i;
   FILE *fp;
   int  nch = 0;
   int  no = 0;
   double   phi[MAXP], psi[MAXP], Eabs[MAXP], Emin, phimin, psimin, dummy;
   double   xphi, xpsi, xEabs;
   char line[MAXLEN+1], *fname;

   /* READ THE REFERENCE MAP */
   Emin = 1E10;
   fname = getstring("Write name on UNIRAS map file: ");
   if ((fp = fileopen(fname, "rt")) != NULL)
   {  printf("UNIRAS map file %s opened for read-only access\n\n", fname);
      printf("\n+---------+---------+---------+---------+\n");
      nch = get_line(fp, line, MAXLEN);
      nch = get_line(fp, line, MAXLEN);
      do 
      {  nch = get_line(fp, line, MAXLEN);
         nch = sscanf(line, "%lf%lf%lf%lf", 
                                      &phi[no], &psi[no], &Eabs[no], &dummy);
         if (Eabs[no] < Emin)
         {  Emin = Eabs[no];
            phimin = phi[no];
            psimin = psi[no];
            putchar('M');
         } 
         else
            putchar('+');
         if (((no+1) % 37) == 0)
            putchar('\n');
         no++;

      } while (nch !=-1);
   }
   fileclose(fname, fp);
   printf("\n+---------+---------+---------+---------+\n");

   printf("\nMINIMUM: PHI= %.2lf  PSI= %.2lf  Emin= %.5lf kcal/mol\n\n",
                           phimin,     psimin,      Emin);

   /* MAIN LOOP */
   do
   {   fname = getstring("Write name on UNIRAS map file: ");
       if ((fp = fileopen(fname, "rt")) != NULL)
       {  printf("UNIRAS map file %s opened for read-only access\n\n", fname);
          printf("\n+---------+---------+---------+---------+\n");
          nch = get_line(fp, line, MAXLEN);
          nch = get_line(fp, line, MAXLEN);
          for(i=0; i<no-1; i++)
          {  nch = get_line(fp, line, MAXLEN);
             nch = sscanf(line, "%lf%lf%lf%lf", &xphi, &xpsi, &xEabs, &dummy);
             if (xEabs < Eabs[i])
             {  Eabs[i] = xEabs;
                phi[i]  = xphi;
                psi[i]  = xpsi;
                if (xEabs < Emin)
                {  Emin = xEabs;
                   phimin = xphi;
                   psimin = xpsi;
                   putchar('M');
                }
                else
                   putchar('+');
             }
             else
                putchar(' ');
             if (((i+1) % 37) == 0)
                putchar('\n'); 
          }
          if (phimin > 180.0)
             phimin -= 360.0;
          if (psimin > 180.0)
             psimin -= 360.0;
          printf("+---------+---------+---------+---------+\n");
          printf("\nMINIMUM: PHI= %.2lf  PSI= %.2lf  Emin= %.5lf kcal/mol\n\n",
                                  phimin,     psimin,      Emin);
          fileclose(fname, fp);
       }
       else
          break;
   } while (fp != NULL);

   printf("\nMINIMUM: PHI= %.2lf  PSI= %.2lf  Emin= %.5lf kcal/mol\n\n",
                           phimin,     psimin,      Emin);
 
   fname = getstring("Write name on UNIRAS target file : ");
   if ((fp = fileopen(fname, "wt")) != NULL)
   {  printf("PHIPSI plot UNIRAS file %s opened for write access\n\n", fname);
      fprintf(fp, "MINIMUM: PHI= %.2lf  PSI= %.2lf  Emin= %.5lf kcal/mol\n",
                                 phimin,     psimin,      Emin);
      fprintf(fp, "     phi       psi            Eabs            Erel\n");
     
      for ( i=0; i<no-1; i++)
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
