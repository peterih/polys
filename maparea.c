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
main()

{  register int i;
   FILE     *fp;
   int      nch = 0;
   int      no = 0;
   double   phi[MAXP], psi[MAXP], Eabs[MAXP], Erel[MAXP];
   double   Emin, Ecut, area, atot, vol3, vol[3];
   double   phimin, psimin;
   char     line[MAXLEN+1], *fname;

   /* READ THE REFERENCE MAP */
   Emin = 1E10;
   fname = getstring("Write name on UNIRAS map file: ");
   if ((fp = fileopen(fname, "rt")) != NULL)
   {  printf("UNIRAS map file %s opened for read-only access\n\n", fname);
      nch = getline(fp, line, MAXLEN);
      nch = getline(fp, line, MAXLEN);
      do 
      {  nch=getline(fp, line, MAXLEN);
         nch=sscanf(line,"%lf%lf%lf%lf",&phi[no],&psi[no],&Eabs[no],&Erel[no]);
         if ((Eabs[no] < Emin) AND (nch !=-1))
         {  Emin = Eabs[no];
            phimin = phi[no];
            psimin = psi[no];
         } 
         no++;
      } while (nch !=-1);
   }
   fileclose(fname, fp);

   printf("\nMINIMUM: PHI= %.2lf  PSI= %.2lf  Emin= %.5lf kcal/mol\n\n",
                           phimin,     psimin,      Emin);

   /* calculate area and 3D volume */
   area = 0.0;
   vol3 = 0.0; 
   vol[0] = 0.0;
   vol[1] = 0.0;
   vol[2] = 0.0;
   Ecut = getreal("Write relative energy cutoff (kcal/mol): ");
   for (i=0; i<no-1; i++)
      if ((phi[i] != 360.00) AND (psi[i] != 360.00))
      {  if ((Eabs[i] - Emin) < Ecut)
         {  area += 100.0; 
            vol3 += Eabs[i]-Emin;
            if ((phi[i] >=   0.0) AND (phi[i] < 120.0))
               vol[0] += Eabs[i]-Emin;
            if ((phi[i] >= 120.0) AND (phi[i] < 240.0))
               vol[1] += Eabs[i]-Emin;
            if ((phi[i] >= 240.0) AND (phi[i] < 360.0))
               vol[2] += Eabs[i]-Emin;
         }
         atot += 100.0;
      }


   printf("         AREA of potential energy surface: %10.3lf\n\n", area);
   printf("   TOTAL AREA of potential energy surface: %10.3lf\n\n", atot);
   printf("OCCUPIED AREA of potential energy surface: %10.3lf\n\n", area/atot);
   printf("    3D VOLUME of potential energy surface: %10.3lf\n\n", vol3*100.0);
   printf("    phi [  0; 120[: %10.3lf\n\n", vol[0]/vol3);
   printf("    phi [120; 240[: %10.3lf\n\n", vol[1]/vol3);
   printf("    phi [240; 360[: %10.3lf\n\n", vol[2]/vol3);

//// 20030117 pih: added 'return' statement in case of correct termination.
   return 0;
} /* End of main */


/* End of file */
