#include <string.h>
#include <math.h>
#include "defs.h"
#include "fileoper.h"

BOOLEAN debug;
/**********************************************************************
	Print welcome header on output device.
**********************************************************************/
main(int argc, char **argv)

{  char    *finame, *foname, str[80];
   int     i, j, idum;
   double  phiarr[1500], psiarr[1500], omearr[1500], earr[1500];  
   double  kT, T, A_0, A_1, Qsum;
   FILE    *fp;

   printf("PROGRAM = %s\n\n", argv[0]);
   if (argc < 3)
   {  printf ("\07");
      printf ("\n\tusage: %s input_file output_file temperature\n\n",argv[0]);
      exit(0);
   }

   if (argc == 4)
      sscanf(argv[3], "%f", &T);
   else
      T  = 298.0;
   printf("TEMPERATURE = %.2f\n\n", T);
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
      {  idum = fscanf(fp, "%lf%lf%lf%lf", &phiarr[i], &psiarr[i], 
                                          &omearr[i], &earr[i]);
         i++;
      } while (idum == 4);

      Qsum = 0.0;
      for (j=0; j<i-1; j++)
      {  A_0 = -earr[j]/kT;
         A_1 = exp(A_0);
         Qsum += A_1;
         earr[j] = Qsum;
         printf("pi[%d] = %.3lf\n", j, earr[j]);
      }
      fileclose(finame,fp);
   }

//// 2003 pih: *foname -> foname . Pointer comparison.
   if (foname != NULL)
   {  if ((fp = fileopen(foname, "wt")) != NULL)
      {  fprintf(fp, "%s\n", str);
         for (j=0; j<i-1; j++)
         {  earr[j] /= Qsum;
            printf("pi[%d] = %.9lf\n", j, earr[j]);
            fprintf(fp, "%7.1lf%7.1lf%7.1lf%20.10lf\n", 
                        phiarr[j], psiarr[j], omearr[j], earr[j]);
         }
         fileclose(foname,fp);
      }
   }

}  /* end of main */
