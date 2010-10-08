/**********************************************************************
  "driver.c"

  This module contains rutines to drive torsional angles through
conformational ("phi-psi") space.

  Written by Soren Balling Engelsen, INRA-93.
**********************************************************************/
#include <stdlib.h>
#include <math.h>
#include "polys.h"
#include "extern.h"


/**********************************************************************
   This rutine will find the network of connected atoms on each
side of the rotatable bond "ai-ae". This is done by virtually deleting
the "ai-ae" bond from the bond table. Then the bond table is passed
once creating a connection level vector CONN[] which finally will
contain *rotlevel for those atoms to be rotated.
   If the specified bond was a part of a cyclic structure then no
action is taken.
**********************************************************************/
int FindConn(int na, int nb, int ai, int ae, int *rotlevel)

{  register int  i; 
   register int  j;
   int      c_ai, c_ae, l_ai, l_ae, level, maxlevel;
   int      bf, bt;

   /* initialize connection level vector */

   for (i=0; i<na; i++)
      CONN[i]=0;

   /* This is the central part of the algorithm which passes through */
   /* the bond table once and generates the connection level vector  */

   maxlevel=0;
   for (j=0; j<nb; j++)
   {  bf = B[j].from;
      bt = B[j].to;
      if (((bf == ai) AND (bt == ae)) OR ((bf == ae) AND (bt == ai)))
      {
      }
      else
      {  maxlevel++;
         if (CONN[bf] != 0) 
         {  level = CONN[bf];
            for (i=0; i<na; i++)
              if (CONN[i] == level)
                  CONN[i] = maxlevel;
         }
         else 
            CONN[bf] = maxlevel;

         if (CONN[bt] != 0) 
         {  level = CONN[bt];
            for (i=0; i<na; i++)
               if (CONN[i]==level)
                  CONN[i] = maxlevel;
         }
         else 
            CONN[bt]=maxlevel;
      }
   }

   l_ai=CONN[ai];
   l_ae=CONN[ae];
   c_ai=0;
   c_ae=0;
   /* If the level of atoms "ai" and "ae" is the same then */
   /* they are connected via another path ==> the bond is  */
   /* a part of a cyclic structure.                        */
   if (l_ai == l_ae) 
   {  for (i=0; i<na; i++)
         CONN[i]=0;
      return(0);
   }

   /* count the number of atoms to be rotated */
   for (i=0; i<na; i++)
   {  if (CONN[i]==l_ai)
         c_ai++;
      if (CONN[i]==l_ae)
         c_ae++;
   }

   if (c_ai <= c_ae)            /* rotate level l_ai */
   {  *rotlevel = l_ai;
      return(c_ai-1);
   }
   else
   {  *rotlevel = l_ae;         /* rotate level l_ae */
      return(-(c_ae-1));
   }
   
} /* End of FindConn */
 

/**********************************************************************
**********************************************************************/
void rotbond(Tors Pi, double torang)

{  int        i, nrot, level;
   double     sign, angle, dangle, tstart; 
   Vector3    vec;
   Matrix     TM;

   /* initializations */
   tstart = cput();
   sign = 1.0;

   /* calculate the present torsional angle */
   angle=torval(A[Pi.from].pos, A[Pi.over1].pos, A[Pi.over2].pos, A[Pi.to].pos);
   dangle = (torang-angle)*DEGtoRAD; 
   if (fabs(dangle) > 0.001)
   {  /* find connection vector */
      nrot = FindConn(M.nat, M.nbond, Pi.over1, Pi.over2, &level);
      /* printf("Rotation about the bond %d - %d include %d atom(s)\n",
                 Pi.over1, Pi.over2, abs(nrot)); */
      if (nrot < 0) 
         sign = -1.0;

      /* calculate the rotation matrix */
      vec.x = A[Pi.over2].pos.x - A[Pi.over1].pos.x;
      vec.y = A[Pi.over2].pos.y - A[Pi.over1].pos.y;
      vec.z = A[Pi.over2].pos.z - A[Pi.over1].pos.z;
      genrot((sign*dangle), A[Pi.over1].pos, vec, TM);

      /* do the rotation */
      for (i=0; i<M.nat; i++)
      {   if (CONN[i] == level)
          {  TransfV(A[i].pos, TM, &vec);
             A[i].pos.x = vec.x;
             A[i].pos.y = vec.y;
             A[i].pos.z = vec.z;
          }
      }
   }
   trot += cput() - tstart;

} /* End of rotbond */
             

/**********************************************************************
   This rutine will drive to torsional angles through a specified
part of conformational space.

option 1:  rigid map (no cartesians saved)
option 2:  rigid map (cartesians saved)
option 3:  rigid map (cartesians saved and full energy evaluation)
option 4:  relaxed map (no cartesians saved)
option 5:  relaxed map (cartesians saved)
option 6:  relaxed map (cartesians saved and full energy evaluation)
**********************************************************************/
void phipsidriv(int phinr, double sphi, double ephi, double incphi,
                int psinr, double spsi, double epsi, double incpsi,
                BOOLEAN silent, int option, char *fname)

{  register int i;  
   int          phicon, psicon;
   double       phiv, cphiv, rphi, Vpot;
   double       psiv, cpsiv, rpsi, hdist;
   double       hfc = 500.0;
   FILE         *fp;

   if (M.nat >= MAXATOMR)
   {  printf("WARNING:\n");
      printf("WARNING: Driver rutine has been setup to molecules less than\n");
      printf("WARNING: %d atoms. No action taken.\n", MAXATOMR);
      printf("WARNING:\n");
   }
   else 
   {  printf("WARNING:\n");
      printf("WARNING: Reference coordinate set overwritten\n");
      printf("WARNING:\n");
   }

   printf("Recording phi-psi Map using the options:\n\n");
   printf("PHI:      %-4s - %-4s - %-4s - %-4s     [%2d]\n",
                     A[P[phinr].from].lab, A[P[phinr].over1].lab,
                     A[P[phinr].over2].lab, A[P[phinr].to].lab, 
                     abs(P[phinr].fold));
   printf("PHI:      %.1lf --[%.1lf]--> %.1lf\n", sphi, incphi, ephi);
   printf("PHI:      Harmonic force constant:   %.2lf kcal/mol\n\n", hfc);
   printf("PSI:      %-4s - %-4s - %-4s - %-4s     [%2d]\n",
                     A[P[psinr].from].lab, A[P[psinr].over1].lab,
                     A[P[psinr].over2].lab, A[P[psinr].to].lab, 
                     abs(P[psinr].fold));
   printf("PSI:      %.1lf --[%.1lf]--> %.1lf\n", spsi, incpsi, epsi);
   printf("PSI:      Harmonic force constraint:   %.2lf kcal/mol\n\n", hfc);

   switch (option)
   { case 1:  printf("Rigid map (no cartesians saved)\n"); break;
     case 2:  printf("Rigid map (cartesians saved)\n"); break;
     case 3:  printf("Rigid map (cartesians saved and full energy evaluation)\n"
);
              break;
     case 4:  printf("Relaxed map (no cartesians saved)\n"); break;
     case 5:  printf("Relaxed map (cartesians saved)\n"); break;
     case 6:  printf("Relaxed map (cartesians saved and full energy evaluation)\
n");
              break;
   }

   phicon = M.nphi+M.ncon;
   psicon = phicon+1;
   M.ncon += 2;

   P[phicon].from = P[phinr].from;
   P[phicon].over1 = P[phinr].over1;
   P[phicon].over2 = P[phinr].over2;
   P[phicon].to = P[phinr].to;
   P[phicon].K1 = hfc;

   P[psicon].from = P[psinr].from;
   P[psicon].over1 = P[psinr].over1;
   P[psicon].over2 = P[psinr].over2;
   P[psicon].to = P[psinr].to;
   P[psicon].K1 = hfc;

   /* save original data in reference coordinate set Ar */
   for (i=0; i<M.nat; i++)
   {  Ar[i].x = A[i].pos.x;
      Ar[i].y = A[i].pos.y;
      Ar[i].z = A[i].pos.z;
   }

   if ((fp = fileopen("PHIPSI.MAP", "wt")) != NULL)

   {  for (phiv = sphi; phiv <= ephi; phiv += incphi)
      {  /* make sure constrain angle is between -180 and 180 */
         cphiv = phiv;
         if (phiv > 180.0) 
            cphiv = phiv - 360.0;
         P[phicon].p1 = DEGtoRAD*cphiv;

         for (psiv = spsi; psiv <= epsi; psiv += incpsi)
         {  /* make sure constrain angle is between -180 and 180 */
            cpsiv = psiv;
            if (psiv > 180.0) 
               cpsiv = psiv - 360.0;
            P[psicon].p1 = DEGtoRAD*cpsiv;
            /* restore original coordinates */
            for (i=0; i<M.nat; i++)
            {  A[i].pos.x = Ar[i].x;
               A[i].pos.y = Ar[i].y;
               A[i].pos.z = Ar[i].z;
            }
            printf("\n--ROT--> PHI=%.2lf  PSI=%.2lf\n", phiv, psiv);
            rotbond(P[phicon], phiv);
            rotbond(P[psicon], psiv);
            if (option < 4)                          /* rigid map */
               Vpot = potener(0);
            else                                     /* relaxed */
               Minimize(0, 10000, 2000, EPSILON, 0.001, TRUE);
            if ((option % 3) == 0) 
               printener(FALSE);
            hdist=dist(A[10].pos,A[32].pos);
/*            fprintf(fp, "--> %10.2lf%10.2lf%20.4lf\n", 
                             phiv, psiv, M.Epot);
*/
            fprintf(fp, "--> %10.2lf%10.2lf%20.4lf%20.10lf\n", 
                             phiv, psiv, M.Epot, hdist);
            rphi = torval(A[P[phicon].from].pos, A[P[phicon].over1].pos,
                          A[P[phicon].over2].pos, A[P[phicon].to].pos);
            rpsi = torval(A[P[psicon].from].pos, A[P[psicon].over1].pos,
                          A[P[psicon].over2].pos, A[P[psicon].to].pos);
            printf("--RES--> PHI=%.2lf  PSI=%.2lf  E=%.4lf\n\n", 
                             rphi,      rpsi,      M.Epot);
            if ((option != 1) AND (option != 4))
               for(i=0; i<M.nat; i++)
                  fprintf(fp, "%15.5lf%15.5lf%15.5lf%5d  %-8s\n", 
                     A[i].pos.x, A[i].pos.y, A[i].pos.z, A[i].ztyp, A[i].lab);

         }
      }
   }
   fileclose("PHIPSI.MAP", fp);

}  /* End of phipsidriv */

/**********************************************************************
   This rutine will drive one torsional angle through a specified
part of conformational space.
**********************************************************************/
void phidriv(int phinr, double sphi, double ephi, double incphi,
             BOOLEAN silent, int option, char *fname)

{  register int i;  
   int          phicon;
   double       phiv, rphi, Vpot;
   double       hfc = 1500.0;
   FILE         *fp;

   if (M.nat >= MAXATOMR)
   {  printf("WARNING:\n");
      printf("WARNING: Driver rutine has been setup to molecules less than\n");
      printf("WARNING: %d atoms. No action taken.\n", MAXATOMR);
      printf("WARNING:\n");
   }
   else 
   {  printf("WARNING:\n");
      printf("WARNING: Reference coordinate set overwritten\n");
      printf("WARNING:\n");
   }

   printf("Recording rotational barrier using the options:\n\n");
   printf("PHI:      %-4s - %-4s - %-4s - %-4s     [%2d]\n",
                     A[P[phinr].from].lab, A[P[phinr].over1].lab,
                     A[P[phinr].over2].lab, A[P[phinr].to].lab, 
                     abs(P[phinr].fold));
   printf("PHI:      %.1lf --[%.1lf]--> %.1lf\n", sphi, incphi, ephi);
   printf("PHI:      Harmonic force constant:   %.2lf kcal/mol\n\n", hfc);
   switch (option)
   { case 1:  printf("Rigid map (no cartesians saved)\n"); break;
     case 2:  printf("Rigid map (cartesians saved)\n"); break;
     case 3:  printf("Rigid map (cartesians saved and full energy evaluation)\n"); 
              break;
     case 4:  printf("Relaxed map (no cartesians saved)\n"); break;
     case 5:  printf("Relaxed map (cartesians saved)\n"); break;
     case 6:  printf("Relaxed map (cartesians saved and full energy evaluation)\n");
              break;
   }
   phicon = M.nphi+M.ncon;
   M.ncon += 1;

   P[phicon].from = P[phinr].from;
   P[phicon].over1 = P[phinr].over1;
   P[phicon].over2 = P[phinr].over2;
   P[phicon].to = P[phinr].to;
   P[phicon].K1 = hfc;

   /* save original data in reference coordinate set Ar */
   for (i=0; i<M.nat; i++)
   {  Ar[i].x = A[i].pos.x;
      Ar[i].y = A[i].pos.y;
      Ar[i].z = A[i].pos.z;
   }

   if ((fp = fileopen("PHI.MAP", "wt")) != NULL)
   {  for (phiv = sphi; phiv <= ephi; phiv += incphi)
      {  P[phicon].p1 = DEGtoRAD*phiv;
         /* restore original coordinates */
         for (i=0; i<M.nat; i++)
         {  A[i].pos.x = Ar[i].x;
            A[i].pos.y = Ar[i].y;
            A[i].pos.z = Ar[i].z;
         }
         printf("\n--ROT--> PHI=%.2lf\n", phiv);
         rotbond(P[phicon], phiv);
         if (option < 4)                          /* rigid map */
            Vpot = potener(0);
         else                                     /* relaxed */
            Minimize(0, 10000, 2000, EPSILON, 0.001, TRUE);
         if ((option % 3) == 0) 
            printener(FALSE);
         fprintf(fp, "--> %10.2lf%20.4lf\n", 
                          phiv, M.Epot);
         rphi = torval(A[P[phicon].from].pos, A[P[phicon].over1].pos,
                       A[P[phicon].over2].pos, A[P[phicon].to].pos);
         printf("--RES--> PHI=%.2lf  E=%.4lf\n\n", 
                          rphi,      M.Epot);
         if ((option != 1) AND (option != 4))
            for(i=0; i<M.nat; i++)
               fprintf(fp, "%15.5lf%15.5lf%15.5lf%5d  %-8s\n", 
                  A[i].pos.x, A[i].pos.y, A[i].pos.z, A[i].ztyp, A[i].lab);

      }
   }
   fileclose("PHI.MAP", fp);

}  /* End of phidriv */


/* End of file */
