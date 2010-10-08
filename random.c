/**********************************************************************
  "random.c"

  This module contains the source code for generating a random
number between 0 and 1. 
  The rutine is implementation specific.

  Written by Soren Balling Engelsen, INRA-93.
**********************************************************************/
#include <stdlib.h>
#include <math.h>
#include "defs.h"

extern BOOLEAN debug;

/**********************************************************************
**********************************************************************/
double get_random(int *idum)

{  static int     iff = 0;
   int            j;
   static double  y, v[98];
   int            dum;

   if ((*idum < 0) OR (iff == 0)) 
   {  iff=1;
      srand(abs(*idum));
      *idum = 1;
      for (j=1; j<=97; j++) 
         dum = rand();
      for (j=1; j<=97; j++) 
         v[j] = (double) rand()/RAND_MAX;
      y = (double) rand()/RAND_MAX;
   }
   j = (int) (1.0 + 97.0*y);
/*   if (debug)
   { printf("\ny = %.3lf\n", y);
     printf("j = %d\n", j);
   }
*/
   if (j > 97 OR j < 1) 
      printf("RANDOM: This cannot happen.");
   y = v[j];
   v[j] = (double) rand()/RAND_MAX;
   return (y);

} /* end of get_random */

/* end of file */
