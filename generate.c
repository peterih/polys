/**********************************************************************
POLYS MODULE <generate.c>

   Written by Soren Balling Engelsen, INRA-93.
**********************************************************************/
#include "polys.h"
#include "extern.h"


/**********************************************************************
   This rutine will generate a table of bonds B[] and check if the 
valencies are saturated.
   Using this table it will generate a table of valence angles (theta) 
T[] and a table of torsional angles (phi) P[] for the polysaccharide
**********************************************************************/
void Generate(int option)

{  double tstart;

   tstart = cput();
   if (M.nbond == 0)
     gen_bondtab();
   gen_angletab();
   gen_tortab();
   gen_nonbond();
   if (option >= 10)
      gen_rotvec();
   tint += cput() - tstart;

}  /* End of Generate */


/* End of file */
