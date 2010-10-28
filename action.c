/**********************************************************************
  "action.c"

  Written by Soren Balling Engelsen, INRA-93, 94.
**********************************************************************/
#include <math.h>
#include <stdio.h>
#include "polys.h"

/**********************************************************************
	Actiontype determines which type of interaction exists between
two atoms a1 and a2. Code for hydrogen bonds not included yet.

return value      interaction type
----------------------------------------
-1                close approach warning
 0                non-bonded
 1                bonded
 2                hydrogen bond and nonbonded

BEWARE! It is not recommended to make any changes here without
making elaborate tests.
**********************************************************************/
int
actiontype(atom a, atom b, double *ab)

#define MAXBDIST   1.81    /* max. bonded distance in topology search */
#define MINBDIST   0.85    /* min. bonded distance in topology search */
#define MAXBHDIST  1.25    /* max. bonded distance to a hydrogen */
#define MINNHDIST  1.15    /* min. bonded distance between non hydrogens */


{  double xd, yd, zd;

   *ab = 9999.0;
   xd = a.pos.x - b.pos.x;
   if (fabs(xd) > MAXBDIST)
      return 0;                      /* non bond */
   yd = a.pos.y - b.pos.y;
   if (fabs(yd) > MAXBDIST)
      return 0;                      /* non bond */
   zd = a.pos.z - b.pos.z;
   if (fabs(zd) > MAXBDIST)
      return 0;                      /* non bond */

   *ab = sqrt(xd*xd + yd*yd + zd*zd);
 
   if (*ab < MAXBHDIST)              /* bond to a hydrogen or WARNING */
   {  /* first general close approach WARNING */
      if (*ab < MINBDIST)
      {  printf("WARNING: atom %-8s[%3d] and atom %-8s[%3d]"
             " is only %5.2f apart.\n\n", a.lab, a.res, b.lab, b.res, *ab);
      }

      /* both are hydrogens -> nonbonded */
      if (((a.ztyp / 100) == 1) AND ((b.ztyp / 100) == 1))
         return 0;

      /* one is a hydrogen */
      if (((a.ztyp / 100) == 1) OR ((b.ztyp / 100) == 1))
         /* only intraresidue covalent bonds to hydrogens */
         if (a.res == b.res)
            return 1;
         else 
            return 0;

      /* close approach WARNING for non-hydrogens */
      if (*ab < MINNHDIST)
      {  if (*ab > MAXBHDIST)
            printf("WARNING: atom %-8s[%3d] and atom %-8s[%3d]"
                " is only %5.2f apart.\n\n", a.lab, a.res, b.lab, b.res, *ab);
         return -1;
      }
   }

   if (*ab <= MAXBDIST) 
   {  if (((a.ztyp / 100) == 1) OR ((b.ztyp / 100) == 1))
         return 2;
      return 1;          /* bond */
   }

   return 0;             /* non bond */

} /* End of actiontype */
