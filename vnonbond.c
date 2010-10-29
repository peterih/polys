/**********************************************************************
 * This module contains the source code for calculation of nonbonded
 * energy as well as its derivatives vrt. cartesians (F). The function
 * includes spherical cutoff (rcut) and minimum imaging in a cubic box.
 * Written by Soren Balling Engelsen, INRA-93.
 */
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "polys.h"
#include "extern.h"

static int nopara;    /* number of nonbonded parameters */
static int nbtype;    /* type of nonbonded potential */
static double dielec; /* the dielectric constant */
static double rcut;   /* the spherical cutoff */

/*
 * This rutine will neutralize overall charge on molecular system.
 * This rutine needs revision to support more than on molecule.
 */
void
neutralize(int nat)
{
   int    i, j, nch;
   double chgsum, correchg;

   chgsum = 0.0;
   for (i=0; i<nat; i++)
      chgsum += A[i].chg;
   printf("Total charge = %.2e\n\n", chgsum);

   if (fabs(chgsum) > TOL)
   {  printf("WARNING: Molecular system is not neutral (%.7lf)\n", chgsum);
      correchg = chgsum/((double) nat);
      for (i=0; i<nat; i++)
         A[i].chg -= correchg;
      printf("WARNING: Atomic (%d) charges added by %.7lf\n", nat, correchg);
      printf("WARNING: Charges specified in input file ");
      printf("has been overwritten\n"); 
   }
}

/**********************************************************************
 * This routine will read nonbonded parameters from the file
 * 'PARAM_DIR/NBPARAM.PF' into the global array NBP[].
 * First it will read some control parameters from the file.
 * ----------------------------------------------------------------
 * nopara : number of parameters to be read (you can have more in the file)
 * nbtype : type of nonbonded potential
 * dielec : dielectric constant
 * rcut   : spherical cutoff distance
 * sc14fv : 14 scaling factor for van der Waals interactions
 * sc14fc : 14 scaling factor for electrostatic interactions
 */
void
read_nbp(int nat, int nbond)
{
   FILE   *fp;
   int    i, typ;
   double Emin, charge, r2, rm, rm6;
   double Aij, Bij, Cij, Eij;
   char   line[MAXLEN];
   char   *pardir = getenv("PARAM_DIR");

   strcpy(line,pardir);
   strcat(line,"/NBPARAM.PF");
   if ((fp = fileopen(line, "rt")) != NULL)
   {  get_line(fp, line, MAXLEN);
      printf("NONBONDED POTENTIAL -> %s\n", line);
      get_line(fp, line, MAXLEN);
      sscanf(line, "%d%d%lf%lf%lf%lf", 
                   &nopara, &nbtype, &dielec, &rcut, &sc14fv, &sc14fc);
      printf("DIELECTRIC CONSTANT     : %.2lf\n", dielec);
      printf("SPHERICAL CUTOFF        : %.2lf\n", rcut);
      printf("1-4 SCALING FACTOR (vdW): %.2lf\n", sc14fv);
      printf("1-4 SCALING FACTOR (Cou): %.2lf\n", sc14fc);
      printf("POTENTIAL TYPE          : ");
 
      switch(nbtype)
      {  case 0: printf("LENNARD-JONES 12-6");
                 break;
         case 1: printf("LENNARD-JONES 12-6-1");
                 break;
         case 2: printf("LENNARD-JONES 12-6-1/charge");
                 break;
         case 3: printf("LENNARD-JONES 12-6-1/dist");
                 break;
         case 4: printf("LENNARD-JONES 9-6-1");
                 break;
         case 5: printf("BUCKINGHAM EXP-6-1");
                 break;
        default: printf("IGNORED");
      }
      printf("\n\n i    ztype       Emin            r2\n");
      printf("=========================================\n");

      for (i=0; i<nopara; i++)
      {  get_line(fp, line, MAXLEN);  
         switch(nbtype)
         { case 0: /* LJ 12-6 Standard (E,r*) */
                   sscanf(line, "%d%lf%lf", &typ, &Emin, &r2);
                   rm = 2.0*r2;
                   rm6 = rm * rm * rm;
                   rm6 *= rm6;
                   NBP[typ].PA = sqrt(4.0*Emin*rm6*rm6);
                   NBP[typ].PC = sqrt(4.0*Emin*rm6);
                   break;
           case 1: /* LJ 12-6-1 Standard (E,r*) */
                   sscanf(line, "%d%lf%lf", &typ, &Emin, &r2);
                   rm = 2.0*r2;
                   rm6 = rm * rm * rm;
                   rm6 *= rm6;
                   NBP[typ].PA = sqrt(4.0*Emin*rm6*rm6);
                   NBP[typ].PC = sqrt(4.0*Emin*rm6);
                   break;
           case 2: /* LJ 12-6-1 Standard (E,r* and charge parameter) */
                   sscanf(line, "%d%lf%lf%lf", &typ, &Emin, &r2, &charge);
                   rm = 2.0*r2;
                   rm6 = rm * rm * rm;
                   rm6 *= rm6;
                   NBP[typ].PA = sqrt(4.0*Emin*rm6*rm6);
                   NBP[typ].PC = sqrt(4.0*Emin*rm6);
                   NBP[typ].PB = charge;
                   break;
           case 3: /* LJ 12-6-1 Nonstandard (E,r*) */
                   sscanf(line, "%d%lf%lf", &typ, &Emin, &r2);
                   NBP[typ].PA = Emin;
                   NBP[typ].PC = r2;
                   break;
           case 4: /* LJ 9-6-1 Standard (E,r*) */
                   sscanf(line, "%d%lf%lf", &typ, &Emin, &r2);
                   NBP[typ].PA = Emin;
                   NBP[typ].PC = r2;
                   break;
           case 5: /* Buckingham exp-6-1 */
                   sscanf(line, "%d%lf%lf", &typ, &Aij, &Bij, &Cij);
                   NBP[typ].PA = Aij;
                   NBP[typ].PB = Bij;
                   NBP[typ].PC = Cij;
                   break;
         }
         printf("%2d: %5d %15.10lf %15.10lf\n", i, typ, Emin, r2);
      }
      /* pair parameters from one atom parameters */
      if (nbtype == 3)
      {
      };
      putchar('\n');
      fileclose("NBPARAM.PF", fp);
      neutralize(nat);
   }
} /* End of read_nbp */

/*
 * This routine gets nonbonded pair parameters in a suitable form from
 * the nonbonded parameter array NBP. This rutine is needed for non-
 * standard Lennard-Jones parameters. NBTYPE == 3 (Amber).
 */
void
get_nbp(int typi, int typj, 
        double *Aij, double *Bij, double *Cij)
{
   register int      i;
   double   r, r3, r6, r12, Emin;
   double   PAi, PAj, PCi, PCj;
   
   PAi = NBP[typi].PA;
   PAj = NBP[typj].PA;
   PCi = NBP[typi].PC;
   PCj = NBP[typj].PC;
   r   = PCi + PCj;
   r3  = r*r*r;
   r6  = r3*r3;
   r12 = r6*r6;
   Emin= sqrt(PAi * PAj);
   *Aij = Emin*r12;
   *Bij = 0.0;
   *Cij = 2.0*Emin*r6;
/*   printf("%3d - %3d   %.9lf   %.9lf\n", typi, typj, *Aij, *Cij); */
}

/* */
double
ANINT(double value)
{
   double  number;

   number = ceil(value);
   if ((value-number) > 0.5) 
      number = number + 1.0; 
   return number;
}

/*
 * Transform from cartesian coordinates to cylindrical coordinates.
 */
void
Cylindrical(Vector3 a1, Vector3 a2, 
            double *r1, double *l1, double *pi,
            double *r2, double *l2, double *pj)
{

}

/*
 * This routine calculates the nonbonded potential. It is prepared
 * for "minimum image convention" and molecular dynamics.
 * derivtyp    0:  no derivatives
 * derivtyp    1:  cartesian derivatives
 * derivtyp    2:  dihedral derivatives
 */
double
Vnonbond(int derivtyp, int nat, double *Vvdw, double *Vc)
{
   register int i;
   register int j;  
   register int k;  
   int      level;
   int      inbex;
   int      i14sc;
   double   tstart;
   double   V, VS;
   double   box;
   double   scvdw, sccou;
   double   Aij, Bij, Cij, Eij;
   double   Qij, Vq, qconst;
   double   boxinv, rcutsq, ex;
   double   lRIJ, RIJsq, irm2, irm6, Vij, Fij, Wij;
   double   li, cri, pi, crj, lj, pj, dvdp, phiv;
   Vector3  FI, RIJ, FIJ;
   atom     AI, AJ;

   /* initialize variables */
   tstart = cput();
   V    = 0.0;
   Vq   = 0.0;
   VS   = 0.0;
   scvdw = 1.0;
   sccou = 1.0;

   /* constants which are needed */
   box    = 999.0;        /* the cubic box length */

   /* calculate useful quantities */
   boxinv = 1.0/box;
   rcutsq = rcut*rcut;
   qconst = 332.091 / dielec;


   /***********************************************************************/
   /* NBTYP = 0   LENNARD JONES 12-6   (standard; no charges; no scaling) */
   /* NBTYP = 1   LENNARD JONES 12-6   (standard; Emin, r2)               */
   /* NBTYP = 2   LENNARD JONES 12-6-1 (standard; Emin, r2; CFF charges)  */
   /* NBTYP = 3   LENNARD JONES 12-6   (nonstandard; pairparam)           */
   /* NBTYP = 4   LENNARD JONES 9-6    (standard; Emin, r2)               */
   /* NBTYP = 5   BUCKINGHAM EXP-6     ()                                 */
   /***********************************************************************/

   /***********************************************************************/
   /* NBTYP = 0   LENNARD JONES 12-6   (standard; no charges; no scaling) */
   /***********************************************************************/
   if (nbtype == 0)
   for (i=0; i < nat-1; i++)
   {  AI = A[i];
      FI  = F[i];
      inbex = INB[i];

      /* Inner loop begins */
      /* This is THE EXPENSIVE LOOP */
      /* Try to avoid nested operations */
      /* Avoid indexing to A[i] */

      for (j=i+1; j < nat; j++)
      {  if ((j != ANB[inbex]) OR (inbex >= INB[i+1]))
         {  RIJ.x = AI.pos.x - A[j].pos.x;
            RIJ.y = AI.pos.y - A[j].pos.y;
            RIJ.z = AI.pos.z - A[j].pos.z;

/*          RIJ.x -= ANINT(RIJ.x*boxinv)*box;
            RIJ.y -= ANINT(RIJ.y*boxinv)*box;
            RIJ.z -= ANINT(RIJ.z*boxinv)*box;                    */

            RIJsq = RIJ.x*RIJ.x + RIJ.y*RIJ.y + RIJ.z*RIJ.z;
            if (RIJsq < rcutsq)
            {  /* Print warning if nb-distance less than 1 A */
               if (RIJsq < 1.0)
                  printf("Atom %s only %.2lf from atom %s\n", 
                          AI.lab, sqrt(RIJsq), A[j].lab);

               /* Get the pair parameters */
               Aij = NBP[AI.ztyp].PA * NBP[A[j].ztyp].PA;
               Bij = NBP[AI.ztyp].PB * NBP[A[j].ztyp].PB;
               Cij = NBP[AI.ztyp].PC * NBP[A[j].ztyp].PC;

               /* calculate the inv-6 factor */
               irm2  = 1.0/RIJsq;
               irm6  = irm2*irm2*irm2;
 
               /* calculate the LJ 12-6 potential */
               Aij  *= irm6;
               Vij   = Aij*irm6 - Cij*irm6;
               V    += scvdw*Vij;

               /* and the derivatives */
               Fij   = Vij + Aij*irm6;
               Fij   *= 6.0 * scvdw / RIJsq;
               FIJ.x = Fij * RIJ.x;
               FIJ.y = Fij * RIJ.y;
               FIJ.z = Fij * RIJ.z;


               FI.x  += FIJ.x;
               FI.y  += FIJ.y;
               FI.z  += FIJ.z;

               F[j].x -= FIJ.x;
               F[j].y -= FIJ.y;
               F[j].z -= FIJ.z;
            }
         }
         else
           inbex++;
      } /* inner loop ends */
       
      F[i] = FI;

   } /* outer loop ends */   

   /***********************************************************************/
   /* NBTYP = 1   LENNARD JONES 12-6   (standard; Emin, r2)               */
   /* NBTYP = 2   LENNARD JONES 12-6-1 (standard; Emin, r2; CFF charges)  */
   /***********************************************************************/
   if ((nbtype == 1) OR (nbtype == 2))
   for (i=0; i < nat-1; i++)
   {  AI = A[i];
      FI  = F[i];
      inbex = INB[i];
      i14sc = I14[i];

      /* Inner loop begins */
      /* This is THE EXPENSIVE LOOP */
      /* Try to avoid nested operations */
      /* Avoid indexing to A[i] */

      for (j=i+1; j < nat; j++)
      {  if ((j != ANB[inbex]) OR (inbex >= INB[i+1]))
         {  RIJ.x = AI.pos.x - A[j].pos.x;
            RIJ.y = AI.pos.y - A[j].pos.y;
            RIJ.z = AI.pos.z - A[j].pos.z;

            RIJsq = RIJ.x*RIJ.x + RIJ.y*RIJ.y + RIJ.z*RIJ.z;
            if (RIJsq < rcutsq)
            {  /* Calculate the nonbonded distance */
               lRIJ = sqrt(RIJsq);

               /* Print warning if nb-distance less than 1 A */
               if (lRIJ < 1.0)
                  printf("Atom %s only %.2lf from atom %s\n", 
                          AI.lab, lRIJ, A[j].lab);

               /* Get the pair parameters */
               Aij = NBP[AI.ztyp].PA * NBP[A[j].ztyp].PA;
               Bij = NBP[AI.ztyp].PB * NBP[A[j].ztyp].PB;
               Cij = NBP[AI.ztyp].PC * NBP[A[j].ztyp].PC;

               /* scaling factor for 1-4 nonbonded interactions  */
               if ((j != A14[i14sc]) OR (i14sc >= I14[i+1]))
               {  sccou = 1.0;
                  scvdw = 1.0;
               }
               else
               {  scvdw = sc14fv;
                  sccou = sc14fc;
                  i14sc++;
               }

               /* calculate the electrostatic potential */
               Qij   = AI.chg*qconst;
               Qij  *= A[j].chg;
               Qij  *= sccou/lRIJ;
               Vq   += Qij;

               /* calculate the 6 factor */
               irm2  = 1.0/RIJsq;
               irm6  = irm2*irm2*irm2;
 
               /* calculate the LJ 12-6 potential */
               Aij  *= irm6;
               Vij   = Aij*irm6 - Cij*irm6;
/*               printf("%d - %d  E = %15.10lf\n", i, j, Vij); */
               V    += scvdw*Vij;

               /* and the derivatives */
               switch (derivtyp)
               {  case 0: break;
                  case 1: Fij   = Vij + Aij*irm6;
                          Fij   *= 6.0 * scvdw / RIJsq;
                          FIJ.x = Fij * RIJ.x;
                          FIJ.y = Fij * RIJ.y;
                          FIJ.z = Fij * RIJ.z;

                          FIJ.x += Qij/RIJsq*RIJ.x;
                          FIJ.y += Qij/RIJsq*RIJ.y;
                          FIJ.z += Qij/RIJsq*RIJ.z;
  
                          FI.x  += FIJ.x;
                          FI.y  += FIJ.y;
                          FI.z  += FIJ.z;
 
                          F[j].x -= FIJ.x;
                          F[j].y -= FIJ.y;
                          F[j].z -= FIJ.z;
                          break;
                  case 2: for (k=0; k<M.nrot; k++)
                          {  /* is this atom going to be rotated */
                             FindConn(M.nat,M.nbond,RT[k].over1,RT[k].over2,&level);
                             if (CONN[j]==level)
                             {  Cylindrical(AI.pos, AJ.pos, &li, &cri, &pi,
                                            &lj, &crj, &pj);
                                phiv = torval(A[RT[k].from].pos, 
                                              A[RT[k].over1].pos,
                                              A[RT[k].over2].pos,
                                              A[RT[k].to].pos);
                                RT[j].p1 = dvdp;
                             }
                          }
                 } /* switch */
             }
         }
         else
           inbex++;
      } /* inner loop ends */
       
      F[i] = FI;

   } /* outer loop ends */   

   /***********************************************************************/
   /* NBTYP = 3   LENNARD JONES 12-6   (nonstandard; pairparam)           */
   /***********************************************************************/
   if (nbtype == 3)
   for (i=0; i < nat-1; i++)
   {  AI = A[i];
      FI  = F[i];
      inbex = INB[i];
      i14sc = I14[i];

      /* Inner loop begins */
      /* This is THE EXPENSIVE LOOP */
      /* Try to avoid nested operations */
      /* Avoid indexing to A[i] */

      for (j=i+1; j < nat; j++)
      {  if ((j != ANB[inbex]) OR (inbex >= INB[i+1]))
         {  RIJ.x = AI.pos.x - A[j].pos.x;
            RIJ.y = AI.pos.y - A[j].pos.y;
            RIJ.z = AI.pos.z - A[j].pos.z;

            RIJsq = RIJ.x*RIJ.x + RIJ.y*RIJ.y + RIJ.z*RIJ.z;
            if (RIJsq < rcutsq)
            {  /* Calculate the nonbonded distance */
               lRIJ = sqrt(RIJsq);

               /* Print warning if nb-distance less than 1 A */
               if (lRIJ < 1.0)
                  printf("Atom %s only %.2lf from atom %s\n", 
                          AI.lab, lRIJ, A[j].lab);

               /* Get the pair parameters */
               get_nbp(AI.ztyp, A[j].ztyp, &Aij, &Bij, &Cij);

               /* scaling factor for 1-4 nonbonded interactions  */
               if ((j != A14[i14sc]) OR (i14sc >= I14[i+1]))
               {  sccou = 1.0;
                  scvdw = 1.0;
               }
               else
               {  scvdw = sc14fv;
                  sccou = sc14fc;
                  i14sc++;
               }

               /* calculate the electrostatic potential */
               Qij   = AI.chg*qconst;
               Qij  *= A[j].chg;
               Qij  *= sccou/lRIJ;
               Vq   += Qij;

               /* calculate the 6 factor */
               irm2  = 1.0/RIJsq;
               irm6  = irm2*irm2*irm2;
 
               /* calculate the LJ 12-6 potential */
               Aij  *= irm6;
               Vij   = Aij*irm6 - Cij*irm6;
               V    += scvdw*Vij;

               /* and the derivatives */
               switch (derivtyp)
               {  case 0: break;
                  case 1: Fij   = Vij + Aij*irm6;
                          Fij   *= 6.0 * scvdw / RIJsq;
                          FIJ.x = Fij * RIJ.x;
                          FIJ.y = Fij * RIJ.y;
                          FIJ.z = Fij * RIJ.z;

                          FIJ.x += Qij/RIJsq*RIJ.x;
                          FIJ.y += Qij/RIJsq*RIJ.y;
                          FIJ.z += Qij/RIJsq*RIJ.z;
  
                          FI.x  += FIJ.x;
                          FI.y  += FIJ.y;
                          FI.z  += FIJ.z;
 
                          F[j].x -= FIJ.x;
                          F[j].y -= FIJ.y;
                          F[j].z -= FIJ.z;
                          break;
                  case 2: for (k=0; k<M.nrot; k++)
                          {  /* is this atom going to be rotated */
                             FindConn(M.nat,M.nbond,RT[k].over1,RT[k].over2,&level);
                             if (CONN[j]==level)
                             {  Cylindrical(AI.pos, AJ.pos, &li, &cri, &pi,
                                            &lj, &crj, &pj);
                                phiv = torval(A[RT[k].from].pos, 
                                              A[RT[k].over1].pos,
                                              A[RT[k].over2].pos,
                                              A[RT[k].to].pos);
                                RT[j].p1 = dvdp;
                             }
                          }
                 } /* switch */
             }
         }
         else
           inbex++;
      } /* inner loop ends */
       
      F[i] = FI;

   } /* outer loop ends */   

   /***********************************************************************/
   /* NBTYP = 4   LENNARD JONES 9-6    (standard; Emin, r2)               */
   /***********************************************************************/
   if (nbtype == 4)
   {
   }

   /***********************************************************************/
   /* NBTYP = 5   BUCKINGHAM EXP-6     ()                                 */
   /***********************************************************************/
   if (nbtype == 5)
   for (i=0; i < nat-1; i++)
   {  AI = A[i];
      FI  = F[i];
      inbex = INB[i];
      i14sc = I14[i];

      /* Inner loop begins */
      /* This is THE EXPENSIVE LOOP */
      /* Try to avoid nested stuff */
      /* Avoid indexing to A[i] */

      for (j=i+1; j < nat; j++)
      {  if ((j != ANB[inbex]) OR (inbex >= INB[i+1]))
         {  RIJ.x = AI.pos.x - A[j].pos.x;
            RIJ.y = AI.pos.y - A[j].pos.y;
            RIJ.z = AI.pos.z - A[j].pos.z;

/*          RIJ.x -= ANINT(RIJ.x*boxinv)*box;
            RIJ.y -= ANINT(RIJ.y*boxinv)*box;
            RIJ.z -= ANINT(RIJ.z*boxinv)*box;                    */

            RIJsq = RIJ.x*RIJ.x + RIJ.y*RIJ.y + RIJ.z*RIJ.z;
            if (RIJsq < rcutsq)
            {  /* calculate the distance */
               lRIJ = sqrt(RIJsq);

               /* print warning if nb-distance less than 1 A */
               if (RIJsq < 1.0)
                  printf("Atom %s only %.2lf from atom %s\n", 
                          AI.lab, sqrt(RIJsq), A[j].lab);

               /* get the pair parameters */
               get_nbp(AI.ztyp, A[j].ztyp, &Aij, &Bij, &Cij);

               /* scaling factor for 1-4 nonbonded interactions  */
               if ((j != A14[i14sc]) OR (i14sc >= I14[i+1]))
               {  sccou = 1.0;
                  scvdw = 1.0;
               }
               else
               {  scvdw = sc14fv;
                  sccou = sc14fc;
                  i14sc++;
               }

               /* calculate the electrostatic potential */
               Qij   = AI.chg*qconst;
               Qij  *= A[j].chg;
               Qij  *= sccou/lRIJ;
               Vq   += Qij;

               /* calculate some EXP-6 factors */
               irm2 = 1.0/RIJsq;
               irm6 = irm2*irm2*irm2;
               ex   = exp(-Bij*lRIJ);
 
               /* calculate the Buckingham EXP-6 potential */
               Aij  *= ex;
               Vij   = Aij - Cij*irm6;
               V    += scvdw*Vij;

               /* and the derivatives */
               Fij   = - Aij*Bij*ex + 6.0*Cij*irm6/lRIJ;
               Fij  *= scvdw/lRIJ;
               FIJ.x = Fij * RIJ.x;
               FIJ.y = Fij * RIJ.y;
               FIJ.z = Fij * RIJ.z;

               FIJ.x += Qij/RIJsq*RIJ.x;
               FIJ.y += Qij/RIJsq*RIJ.y;
               FIJ.z += Qij/RIJsq*RIJ.z;

               FI.x  += FIJ.x;
               FI.y  += FIJ.y;
               FI.z  += FIJ.z;

               F[j].x -= FIJ.x;
               F[j].y -= FIJ.y;
               F[j].z -= FIJ.z;
            }
         }
         else
           inbex++;
      } /* inner loop ends */
       
      F[i] = FI;

   } /* outer loop ends */   

   /* calculate shifted potential  */
/*   sr2 = sigsq / rcutsq;
   sr6 = sr2*sr2*sr2;
   sr12= sr6*sr6;
   Vij = sr12 - sr6;
   VS  = V - ((double) ncut) * Vij;  */

   /* return values */
   *Vvdw = V;
   *Vc   = Vq;
   tnb += cput() - tstart;
   return VS;
}  /* End of Vnonbond */
