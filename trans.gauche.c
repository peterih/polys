#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "defs.h"
#include "fileoper.h"

BOOLEAN debug;
/**********************************************************************
	Print welcome header on output device.
**********************************************************************/
main(int argc, char **argv)

{  char    *finame, *foname, str[80], dum[4];
   int     i, j, k, idum;
   double  phiarr[1500], psiarr[1500], omearr[1500], earr[1500];  
   double  kT, T, A_0, A_1, Qsum, piall, pitrans, pigauche, Emin=1E6;
   FILE    *fp;

   printf("PROGRAM = %s\n\n", argv[0]);

   if (argc < 3)
   {  printf("\n%s\n", "Usage:");
      printf("\nUsage: %s <input> <output> <temp>\n");
      printf("\nUsage:\n");
      exit(0);
   }

   /* TEMPERATURE is argv[3] */
   idum = sscanf(argv[3], "%lf", &T);
   if (idum < 0)
      T = 298.0;
   printf("TEMPERATURE = %lf\n\n", T);
   kT = (1.38062E-23 * T * 6.0221E23)/4.1855E3;

   finame = argv[1];
   printf("INPUT FILE NAME = %s\n\n", finame);
   foname = argv[2];
   printf("OUTPUT FILE NAME = %s\n\n", foname);

   if ((fp = fileopen(finame, "rt")) != NULL)
   {  idum = get_line(fp, str, 80);
      printf("INFO: %s\n\n", str);

      i = 0;
      do
      {  idum = fscanf(fp, "%s%lf%lf", dum, &phiarr[i], &earr[i]);
         if (earr[i] < Emin)
            Emin = earr[i];
         i++;
      } while (idum == 3);

      fileclose(finame,fp);
   }
   /* store i */

   /* let get relative energies */
   for (k=0; k<i-1; k++)
      earr[k] -= Emin;
   

   Qsum = 0.0;
   for (j=0; j<i-1; j++)
   {  A_0 = -earr[j]/kT;
      A_1 = exp(A_0);
      Qsum += A_1;
      earr[j] = A_1;
   }

   /* 2003 pih: *foname -> foname . Pointer comparison */
   if (foname != NULL)
   {  if ((fp = fileopen(foname, "wt")) != NULL)
      {  fprintf(fp, "%s\n", str);
         for (j=0; j<i-1; j++)
         {  earr[j] /= Qsum;
            printf("pi[%5.1lf] = %.5lf\n", phiarr[j], earr[j]);
            fprintf(fp,"%7.1lf%20.10lf\n", phiarr[j], earr[j]);
         }
         
         piall = 0.0;
         pitrans = 0.0;
         pigauche = 0.0;
         for (k=0; k<i-1; k++)
         {   piall += earr[k];
             if (fabs(phiarr[k]) > 120.0)
                pitrans += earr[k];
             else 
                pigauche += earr[k];
         }
         fprintf(fp, "\n\nSUM(pi)         = %.5lf\n", piall);
         fprintf(fp, "trans fraction  = %.1lf\n", pitrans*100.0);
         fprintf(fp, "gauche fraction = %.1lf\n", pigauche*100.0);
         printf("\n\nSUM(pi)         = %.5lf\n", piall);
         printf("trans fraction  = %.1lf\n", pitrans*100.0);
         printf("gauche fraction = %.1lf\n", pigauche*100.0);
         fileclose(foname,fp);
      }
   }
}  /* end of main */
