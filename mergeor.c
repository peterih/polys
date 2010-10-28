/**********************************************************************
  "mergeor.c"

  This program takes a number of optical rotation map files and merge
then into a single "adiabatic" map file.

  Written by Soren Balling Engelsen, INRA-93.
**********************************************************************/
#include <math.h>
#include <string.h>
#include "defs.h"
#include "sinput.h"
#include "fileoper.h"
 
/**********************************************************************
        Reads optical rotation *.or files
**********************************************************************/
BOOLEAN  debug;

#define    MAXP   1500

//// 20030117 pih: removed 'void' declaration of main.
main()

{  register int i;
   FILE *fp;
   int  nch = 0;
   int  no = 0;
   double   phi[MAXP], psi[MAXP], Eabs[MAXP], Rot[MAXP];
   double   Emin, GEmin, phimin, psimin, rotmin, dummy;
   double   xphi, xpsi, xEabs, xRot;
   char     line[MAXLEN+1], *fname;

   /**************************************************************
       READ THE REFERENCE MAP 
   **************************************************************/
   Emin = 10000.0;
   GEmin = getreal("Write global energy minimum: ");
   fname = getstring("Write name on Optical Rotation map file: ");
   if ((fp = fileopen(fname, "rt")) != NULL)
   {  printf("OR map file %s opened for read-only access\n\n", fname);
      printf("\n+---------+---------+---------+---------+\n");
      do 
      {  nch = get_line(fp, line, MAXLEN);
         nch = sscanf(line, "%lf%lf%lf%lf%lf%lf", 
              &phi[no], &psi[no], &Eabs[no], &dummy, &dummy, &Rot[no]);

         /* get the "best" values */
         if ((Eabs[no] < Emin) AND (nch == 6))
         {  Emin = Eabs[no];
            printf("Eminr = %.3lf no = %d nch = %d \n", Emin, no, nch);
            phimin = phi[no];
            psimin = psi[no];
            rotmin = Rot[no];
            putchar('M');
         } 
         else
            putchar('+');
         if (((no+1) % 37) == 0)
            putchar('\n');
         no++;

      } while (nch == 6);
   }
   fileclose(fname, fp);
   printf("\n+---------+---------+---------+---------+\n");

   printf("\nMINIMUM: PHI= %.2lf  PSI= %.2lf  Emin= %.5lf kcal/mol OR= %.2lf\n\n",
                           phimin,     psimin,      Emin, rotmin);

   /* MAIN LOOP */
   do
   {   fname = getstring("Write name on Optical Rotation map file: ");
       if ((fp = fileopen(fname, "rt")) != NULL)
       {  printf("OR map file %s opened for read-only access\n\n", fname);
          printf("\n+---------+---------+---------+---------+\n");
          for(i=0; i<no-1; i++)
          {  nch = get_line(fp, line, MAXLEN);
             nch = sscanf(line, "%lf%lf%lf%lf%lf%lf", 
                   &xphi, &xpsi, &xEabs, &dummy, &dummy, &xRot);

             if (xEabs < Eabs[i])
             {  Eabs[i] = xEabs;
                phi[i]  = xphi;
                psi[i]  = xpsi;
                Rot[i]  = xRot; 

                if (xEabs < Emin)
                {  Emin = xEabs;
                   phimin = xphi;
                   psimin = xpsi;
                   rotmin = xRot;
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
          printf("+---------+---------+---------+---------+\n");
          printf("\nMINIMUM: PHI= %.2lf  PSI= %.2lf  Emin= %.5lf kcal/mol OR= %.2lf\n\n",
                                  phimin,     psimin,      Emin, rotmin);
          fileclose(fname, fp);
       }
       else
          break;
   } while (fp != NULL);

   printf("\nMINIMUM: PHI= %.2lf  PSI= %.2lf  Emin= %.5lf kcal/mol OR= %.2lf\n\n",
                           phimin,     psimin,      Emin, rotmin);
 
   fname = getstring("Write name on Optical Rotation target file : ");
   if ((fp = fileopen(fname, "wt")) != NULL)
   {  printf("Merged OR file %s opened for write access\n\n", fname);
      fprintf(fp, "MINIMUM: PHI= %.2lf  PSI= %.2lf  Emin= %.5lf kcal/mol OR= %.2lf\n",
                                 phimin,     psimin,      Emin, rotmin);
      fprintf(fp, "     phi       psi            Erel            OR\n");
     
      for ( i=0; i<no-1; i++)
      {  if ((Eabs[i] - GEmin) > 20.0)
         {  Eabs[i] = Emin + 20.0;
            Rot[i] = -200.0;
         }
         if (phi[i] < 0.0)
            phi[i] += 360.0;
         if (psi[i] < 0.0)
            psi[i] += 360.0;
         fprintf(fp, "%10.1lf%10.1lf%15.2lf%15.2lf\n",
                      phi[i], psi[i], Eabs[i] - Emin, Rot[i]);
      }
   }
   fileclose(fname, fp);

//// 20030117 pih: added 'return' statement in case of correct termination.
   return 0;
} /* End of main */


/* End of file */
