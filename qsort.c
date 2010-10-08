/**********************************************************************
  "Qsort.c"

  Written by Soren Balling Engelsen, INRA-93.
**********************************************************************/
#include "defs.h"

void
Qsort(int low, int high, int ia[]) 
{  
   int   i, j, midt, save;

   if (low < high) {
      i = low;
      j = high;
      midt = ia[(i+j) / 2];
      do
      {  while (ia[i] < midt) i++;
         while (ia[j] > midt) j--;
         if (i <= j) 
         {  SWAP(save, ia[i], ia[j]);
            i++;
            j--;
         }
      } while (j >= i);
      Qsort(low, j, ia);
      Qsort(i, high, ia);
   }   
} /* End of Qsort */
