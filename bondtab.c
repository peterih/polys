/**********************************************************************
  "bondtab.c"

  Written by Soren Balling Engelsen, INRA-93, 94.
**********************************************************************/
#include <stdio.h>
#include <math.h>
#include "polys.h"
#include "extern.h"


/**********************************************************************
	Returns the valency of atom A
**********************************************************************/
int Valency(atom A)
 
{  return (A.ztyp % 10);

}  /* End of Valency */

/**********************************************************************
	Calculates the distance between two points in cartesian space.
**********************************************************************/
double dist(Vector3 a, Vector3 b)

{  double dx, dy, dz;

   dx = a.x - b.x;
   dy = a.y - b.y;
   dz = a.z - b.z;
   return(sqrt(dx*dx + dy*dy + dz*dz));

} /* End of dist */

/**********************************************************************
**********************************************************************/
void add_bond(int i, int j)

{  B[M.nbond].from=i, B[M.nbond].to=j;
   M.nbond++;

} /* End of add_bond */

/**********************************************************************
	Generates bondtabel B[] and counts all bonds, nonbonds, close
approach warnings and hydrogen bonds. The routine also checks that all
valencies are filled.
**********************************************************************/
void gen_bondtab()

{  register int i;
   register int j;
   int      hb = 0;
   int      nnb = 0;
   int      nw = 0;
   int      nosvw = 0;    /* number of oversaturated valence warnings */
   int      nusvw = 0;    /* number of undersaturated valence warnings */
   int      act, close;
   int      valence, vali, valj, valtab[MAXATOM];
   double   dij, minnbdist;

   printf("Generating bond tabel for %s\n\n", M.id);
   M.nbond=0;                      /* reset number of bonds */
   for ( i=0; i<M.nat; i++)
      valtab[i] = 0;
 
   for ( i=0; i<M.nat; i++)
      for (j=i+1; j<M.nat; j++)
      {  act = actiontype(A[i], A[j], &dij);
         switch(act)
         { case 0:  nnb++;
                    break;
           case 1:  B[M.nbond].from=i, B[M.nbond].to=j;
                    M.nbond++;
                    valtab[i]++, valtab[j]++;
                    break;
           case -1: nw++;
                    break;
           case 2:  nnb++;
                    hb++;
                    break;
         }
      }
   
   for ( i=0; i<M.nat; i++)
   {  valence = (A[i].ztyp % 10);
      if ((valtab[i] < valence) AND (valence != 9))
      {  printf("WARNING: Unfilled valency for atom --> %-8s[%3d] (%d)\n\n",
                   A[i].lab, A[i].res, valtab[i]);
         nusvw++;
      }
      if (valtab[i] > valence)
      { printf("WARNING: Oversaturated valency for atom --> %-8s[%3d] (%d)\n\n",
                   A[i].lab, A[i].res, valtab[i]);
         nosvw++;
      }
   }

   printf("  Number of atoms                          : %7d\n", M.nat);
   printf("  Number of bonded interactions            : %7d\n", M.nbond);
   printf("  Number of non-bonded interactions        : %7d\n", nnb);
   printf("  Number of hydrogen bonds                 : %7d\n", hb);
   printf("  Number of close approach warnings        : %7d\n", nw);
   printf("  Number of oversaturated valency warnings : %7d\n", nosvw);
   printf("  Number of undersaturated valency warnings: %7d\n\n", nusvw);

   /* create missing bonds */
   if ((nusvw != 0) AND (nosvw == 0))
   {  printf("BEWARE: POLYS will try to create missing bonds\n\n");
      for ( i=0; i<M.nat; i++)
      {  vali = (A[i].ztyp % 10);
         if ((valtab[i] < vali) AND (vali != 9))
         {  printf("ATOM %-8s [%d]: has %d unsaturated valency\n",
                    A[i].lab, A[i].res, vali-valtab[i]);
            minnbdist = 9999.0; 
            for ( j=i+1; j<M.nat; j++)
            {  if (actiontype(A[i], A[j], &dij) != 1)  /* not-bonded */
               {  if ((A[i].res == A[j].res) AND (dij < minnbdist))
                  {  close = j;
                     minnbdist = dij;
                  }
               }
            } 
            valj = (A[j].ztyp % 10);
            if ((valtab[j] < valj) AND (minnbdist < 2.0))
            {  B[M.nbond].from=i, B[M.nbond].to=close;
               M.nbond++;
               valtab[i]++, valtab[j]++;
               printf("Abnormal bond length found:\n");
               printf("%-8s[%3d]  <-- %.2lf -->  %-8s[%3d]\n\n",
                   A[i].lab, A[i].res, minnbdist, A[close].lab, A[close].res);
            }
            else
            {  printf("WARNING\n");
               printf("WARNING: vacant valency could not be filled\n");
               A[i].ztyp--;
               printf("WARNING: Atom type for %s changed to %d\n",
                        A[i].lab, A[i].ztyp);
               printf("WARNING\n");
            }
         }
      }
   }
            
}  /* End of gen_bondtab */

/* End of file */
