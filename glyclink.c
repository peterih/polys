/**********************************************************************
POLYS MODULE <glyclink.c>

  This module contains the source code for generating the glycosidic 
linkage between two segments of the progressing polysaccharide chain.

  Written by Soren Balling Engelsen, INRA-93.
  Last revision: S.B.E. 1997
**********************************************************************/
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>
#include "polys.h"
#include "extern.h"


/**********************************************************************
    This routine finds and return the atom index by comparing the 
label "lab" with A[].lab. The rutine is thus Label specific. The rutine
search backwards since the polymer is most likely to be prolonged from
the previous residue.
**********************************************************************/
int AtomIndex(int resnr, char *lab) 

{  int i;

   for (i=M.nat-1; i>=0; i--)
      if (A[i].res == resnr)
         if (strcmp(lab, A[i].lab) == 0)
            return i;
   printf("FindAIndex: Unable to find atom %8s[%3d]\n", lab, resnr);
   printf("FindAIndex: Last Atom %8s[%3d]\n", A[M.nat-1].lab, A[M.nat-1].res);
   exit(0);
 
} /* End of AtomIndex */


/**********************************************************************
    This routine deletes the atom with the index 'no' from the atom
array A[]. All atom indicies above 'no' is shiftet minus one.
**********************************************************************/
void DelAtom(int no) 

{  register int i;

   if (debug)
      printf("DELATOM: Deleting atom #%d\n", no);
   for (i=no+1; i<M.nat; i++)
      A[i-1] = A[i];
   M.nat--;
 
} /* End of DelAtom */

/**********************************************************************
    This routine deletes the Hydrolysis DA array at the end of the day.
The sorting and reverse deletion of the DA array is necessary in order 
to memorize the correct numbers. The variable nwa is the number of
"water" atoms to be deleted (usually the whole array - da).
**********************************************************************/
void Hydrolysis(int nwa) 

{  register int i,j;

   if (nwa>0)
   {  Qsort(da-nwa,da-1,DA);
      for (i=da-1; i>=da-nwa; i--)
      {  /* printf("Deleting atom: %d - %s\n", DA[i], A[DA[i]].lab); */
         DelAtom(DA[i]);
         /* remember to shift other atomnumbers higher than DA[i] */
         for (j=0; j<da-nwa; j++)
            if (DA[j] > DA[i])
               DA[j]--;
      }
      da-=nwa;
   }
 
} /* End of Hydrolysis */


/**********************************************************************
**********************************************************************/
void getconval(char *fname, double *phip, double *psip, double *omep)

{  FILE    *fp;
   char    line[MAXLEN+1];
   int     idum = (-1);
   double  pr, pi, oldpi, phi, psi, ome;


   *phip = 999.9;
   *psip = 999.9;
   *omep = 999.9;

   DelQuotes(fname);
   pr = get_random(&seed); 
   if (debug) printf("\nGETCONVAL: #RAND = %.6lf\n", pr);

   if ((fp = fileopen(fname, "rt")) != NULL)
   {  idum = getline(fp, line, MAXLEN);

      oldpi=0.0;
      for (;(idum != -1);)
      {  idum = getline(fp, line, MAXLEN);
         idum = sscanf(line, "%lf%lf%lf%lf", &phi, &psi, &ome, &pi);
         if INRANGE(oldpi, pi, pr)
         {  *phip = phi;
            *psip = psi;
            *omep = ome;
            break;
         }
         oldpi = pi;
      }
      fileclose(fname, fp);
   }
   else
   {  printf("\n");
      printf("Could not find linkage information in file: %s\n", fname);
      printf("\n");
      exit(0);
   }

} /* End of getconval */

/**********************************************************************
**********************************************************************/
void TransfSeg(int res1, int res2, Matrix TM)

{  int      i;
   Vector3 av;

   for ( i=0; ((i<M.nat) AND (A[i].res<=res2)); i++)
      if (A[i].res >= res1)
      {   av = A[i].pos;
          TransfV(av, TM, &A[i].pos);
      }

} /* End of Transform */

/**********************************************************************
   This rutine connects two saccharide segments by a glycosidic
linkage. Both segments must have been read into the atom array before 
this rutine is called. Hydrolysis must be done afterwards due to the
generality of the rutine, but indicies of all atoms affected is stored 
in a one-dimensional vector DA[] for later deletion.
----------------------------------------------------------
red1    : first residue on the reducing segment
redres  : reducing residue to connect
red2    : last residue on the reducing segment
nred1   : first residue on the non-reducing segment
nredres : non-reducing residue to connect
nred1   : last residue on the non-reducing segment
redmax  :  
redcon  : no. of connection hydroxyl group on reducing residue
nredcon : no. of connection hydroxyl group on non-reducing residue
phi     : value of the desired phi angle
psi     : value of the desired psi angle
omega   : value of the desired ome angle
----------------------------------------------------------
**********************************************************************/
void GlycLink(int red1, int redres, int red2, 
              int nred1, int nredres, int nred2,
              int redmax, int redcon, int nredcon, 
              double phi, double psi, double omega)

{  int     i, HOr, On, HOn, Ar, COr, Or, COn, An, An2, HR, HS;
   double  tau, ang_b, ang_a;
   Vector3 vec, vecHR, vecHS;
   Matrix  TM;

   if (debug)
   {  printf("GLYCLINK: red1  = %d  redres  = %d  red2  = %d  ", 
                                red1,         redres,     red2);
      if (Res[redres].ringtype==2)
        printf("cetofuranose\n");
      if (Res[redres].ringtype==5)
        printf("aldofuranose\n") ;
      if (Res[redres].ringtype==6)
        printf("pyranose\n") ;
      if (Res[redres].ringtype==7)
        printf("heptose\n");
      printf("GLYCLINK: nred1 = %d  nredres = %d  nred2 = %d  ", 
                                nred1,        nredres,    nred2);
      if (Res[nredres].ringtype==2)
        printf("cetofuranose\n");
      if (Res[nredres].ringtype==5)
        printf("aldofuranose\n");
      if (Res[nredres].ringtype==6)
        printf("pyranose\n");
      if (Res[nredres].ringtype==7)
        printf("heptose\n");
   }

   /* 
   find atomic indicies for the atoms involved in the linkage
      DELETE:        HOr, On, HOn
      PHI:           Ar - COr - Or - COn
      PSI:           COr - Or - COn - An
      OMEGA:         Or - COn - An - An2
      TAU:           COr - Or - COn
   */

   /* reducing end */
   switch(redcon)
   {  case 1: if (Res[redres].ringtype==6)         /* pyranose */
              {  Ar  = AtomIndex(redres, "O_5");
                 COr = AtomIndex(redres, "C_1");
                 Or  = AtomIndex(redres, "O_1");
                 HOr = AtomIndex(redres, "HO_1");
              }
              else if (Res[redres].ringtype==5)    /* furanose */
              {  Ar  = AtomIndex(redres, "O_4");
                 COr = AtomIndex(redres, "C_1");
                 Or  = AtomIndex(redres, "O_1");
                 HOr = AtomIndex(redres, "HO_1");
              }
              else
              {
                 printf("POLYS: Unknown linkage type\n");
              }
              break;
      case 2: if (Res[redres].ringtype==7)
              {  Ar  = AtomIndex(redres, "O_6");      /* heptose */
                 COr = AtomIndex(redres, "C_2");
                 Or  = AtomIndex(redres, "O_2");
                 HOr = AtomIndex(redres, "HO_2");
              }
              else
              {
                 printf("POLYS: Unknown linkage type\n");
              }
              break;
   default:printf("\nGlycLink: Glucosidic linkage between residue %d\n",redres);
              printf("GlycLink: and residue %d NOT possible\n",nredres);
              printf("GlycLink: Connection to reducing segment: %d\n\n\n",
                         redcon);
              exit(0);
   }
   /* nonreducing end */
   switch(nredcon)
   {  case 1: HOn = AtomIndex(nredres, "HO_1");
              On  = AtomIndex(nredres, "O_1");
              COn = AtomIndex(nredres, "C_1");
              An  = AtomIndex(nredres, "O_5");
              tau = 116.0;
              printf("BEWARE non-reducing linkage. Please check\n");
              break;
      case 2: HOn = AtomIndex(nredres, "HO_2");
              On  = AtomIndex(nredres, "O_2");
              COn = AtomIndex(nredres, "C_2");
              An  = AtomIndex(nredres, "C_3");
              tau = 116.0;
              break;
      case 3: HOn = AtomIndex(nredres, "HO_3");
              On  = AtomIndex(nredres, "O_3");
              COn = AtomIndex(nredres, "C_3");
              An  = AtomIndex(nredres, "C_4");
              tau = 116.0;
              break;
      case 4: HOn = AtomIndex(nredres, "HO_4");
              On  = AtomIndex(nredres, "O_4");
              COn = AtomIndex(nredres, "C_4");
              An  = AtomIndex(nredres, "C_5");
              tau = 116.5;
              break;
      case 5: if (Res[nredres].ringtype==5)         /* aldofuranose only */
              {  HOn = AtomIndex(nredres, "HO_5"); 
                 On  = AtomIndex(nredres, "O_5");
                 COn = AtomIndex(nredres, "C_5");
                 HR  = AtomIndex(nredres, "H_51_R");
                 HS  = AtomIndex(nredres, "H_52_S");
                 An  = AtomIndex(nredres, "C_4");
                 An2 = AtomIndex(nredres, "O_4");
                 tau = 114.5;
              }
              else if (Res[nredres].ringtype==7)   /* heptose only */
              {  HOn = AtomIndex(nredres, "HO_5");
                 On  = AtomIndex(nredres, "O_5");
                 COn = AtomIndex(nredres, "C_5");
                 An  = AtomIndex(nredres, "C_6");
                 tau = 116.5;
              }
              else
              {
                 printf("POLYS: Unknown linkage type\n");
              }
              break;
      case 6: HOn = AtomIndex(nredres, "HO_6");
              On  = AtomIndex(nredres, "O_6");
              COn = AtomIndex(nredres, "C_6");
              HR  = AtomIndex(nredres, "H_61_R");
              HS  = AtomIndex(nredres, "H_62_S");
              An  = AtomIndex(nredres, "C_5");
              An2 = AtomIndex(nredres, "O_5");
              tau = 114.5;
              break;
      case 8: HOn = AtomIndex(nredres, "HO_8");
              On  = AtomIndex(nredres, "O_8");
              COn = AtomIndex(nredres, "C_8");
              An  = AtomIndex(nredres, "C_7");
              An2 = AtomIndex(nredres, "O_7");
              tau = 118.0;                      /* NeuAc a2-8 NeuAc */
              break;
     default:printf("GLYCLINK: Glucosidic linkage between residue %d\n",redres);
               printf("GLYCLINK: and residue %d NOT possible\n", nredres);
               printf("GLYCLINK: Connection to non-reducing segment: %d\n\n\n",
                         nredcon);
               exit(0);
   }

   /* Orient segment with the reducing connection */
   Orient(red1, red2, Or, HOr, COr);
   /* Orient(1, redmax, Or, HOr, COr);  */

   /* Orient segment with the non-reducing connection */
   Orient(nred1, nred2, HOn, On, COn);

   /* translate segment with nonreducing segment in order to
      have On overlapping Or in origo                  */ 
   tran3(A[On].pos.x, A[On].pos.y, A[On].pos.z, TM);
   TransfSeg(nred1,nred2,TM);

   /*************************************************************/ 
   /* rotate nonreducing part to adjust COr-Or-COn angle to tau */
   /*************************************************************/ 
   ang_b = angval(A[COr].pos, A[Or].pos, A[COn].pos);
   rot3(2, (tau - ang_b)*DEGtoRAD, TM);
   TransfSeg(nred1, nred2, TM); 
   ang_a = angval(A[COr].pos, A[Or].pos, A[COn].pos);
   Res[redres].theta = ang_a;
   if (debug) 
     printf("GLYCLINK: TAU: %8.2lf   -->   %8.2lf\n", ang_b, ang_a);
   
   /*************************************************************/ 
   /* adjust phi-torsional angle -  PHI:   Ap - COp - Op - COm  */
   /*************************************************************/ 
   ang_b = torval(A[Ar].pos, A[COr].pos, A[Or].pos, A[COn].pos);
   vec.x = A[COr].pos.x - A[Or].pos.x;
   vec.y = A[COr].pos.y - A[Or].pos.y;
   vec.z = A[COr].pos.z - A[Or].pos.z;
   genrot((phi - ang_b)*DEGtoRAD, A[Or].pos, vec, TM);
   TransfSeg(nred1, nred2, TM);
   ang_a = torval(A[Ar].pos, A[COr].pos, A[Or].pos, A[COn].pos);
   if (debug) 
     printf("GLYCLINK: PHI: %8.2lf   -->   %8.2lf\n", ang_b, ang_a);
    
   /*************************************************************/ 
   /* adjust psi-torsional angle -  PSI:  COr - Or - COn - An   */
   /*************************************************************/ 
   ang_b = torval(A[COr].pos, A[Or].pos, A[COn].pos, A[An].pos);
   vec.x = A[Or].pos.x - A[COn].pos.x;
   vec.y = A[Or].pos.y - A[COn].pos.y;
   vec.z = A[Or].pos.z - A[COn].pos.z;
   genrot((psi - ang_b)*DEGtoRAD, A[COn].pos, vec, TM);
   TransfSeg(nred1, nred2, TM);
   ang_a = torval(A[COr].pos, A[Or].pos, A[COn].pos, A[An].pos);
   if (debug) 
     printf("GLYCLINK: PSI: %8.2lf   -->   %8.2lf\n", ang_b, ang_a);
    
   /*************************************************************/ 
   /* if 1-5, 1-6 or 2-8 linkage adjust omega-torsional angle   */
   /*       OMEGA:  Or - COn - An - An2                         */
   /*************************************************************/ 
   if ((nredcon >= 5) AND (Res[nredres].ringtype<7))
   {  ang_b = torval(A[Or].pos, A[COn].pos, A[An].pos, A[An2].pos);
      vec.x = A[COn].pos.x - A[An].pos.x;
      vec.y = A[COn].pos.y - A[An].pos.y;
      vec.z = A[COn].pos.z - A[An].pos.z;
      genrot((omega - ang_b)*DEGtoRAD, A[An].pos, vec, TM);
      /* save positions of methylene hydrogens */
      vecHR = A[HR].pos;
      vecHS = A[HS].pos;
      TransfSeg(nred1, nred2, TM);
      ang_a = torval(A[Or].pos, A[COn].pos, A[An].pos, A[An2].pos);
      if (debug) 
        printf("GLYCLINK: OME: %8.2lf   -->   %8.2lf\n", ang_b, ang_a);
      /* restore positions of methylene hydrogens */
      A[HR].pos=vecHR;
      A[HS].pos=vecHS;
   }

   /* Change atom type of Or to 'ether' oxygen   */
   /* Important when using PEF***                */
   A[Or].ztyp = 812;

   /* finally update the hydrolysis array */
   DA[da++]=HOn;
   DA[da++]=On;
   DA[da++]=HOr;
/* for (i=0; i<da; i++)
        printf("DA[%d] = %d\n",i,DA[i]);  */

}  /* End of GlycLink */

/**********************************************************************
   This rutine repeats a segment one time by a simple copy procedure.
   The segment to be repeated must have been read into the atom array
   before the rutine is called. 
----------------------------------------------------------
rep     : the number of repeats
bres    : the first residue of the repeat segment
eres    : the last residue of the repeat segment
----------------------------------------------------------
**********************************************************************/
void RepeatStruc(int rep, int bres, int eres)

{  int  i, j, newat, offres;
   BOOLEAN   first;

   printf("Repeating unit %d of the chain between ", rep);
   printf("residue %d and residue %d\n", bres, eres);
   printf("By copying all residues from %d to %d\n", bres, eres);

   Qsort(0,da-1,DA);
   first = TRUE;
   offres=0;
   newat=0;
   j=0;
   printf("Nnat = %d\n", M.nat);
   printf("Nnres = %d\n", M.nres);
   for (i=0; ((i<M.nat) AND (A[i].res <= eres)); i++)
      if (A[i].res >= bres)
      {  if (first)
         {  M.nres++;
            first = FALSE;
         }
         else
            if (A[i].res != A[i-1].res)
               M.nres++;
         A[M.nat+newat] = A[i];
         A[M.nat+newat].res = M.nres - 1;
/*         if ((i==DA[j]) AND (da > 0))
         {  DA[da++]=M.nat+newat;
            printf("DA[%d] = %d\n",da-1,DA[da-1]);
            j++;
         }                 */
         newat++;
      }
/*   for (i=0; i<da; i++)
      printf("DA[%d] = %d\n",i,DA[i]);   */
       
   M.nat += newat;
   printf("Nnat = %d\n", M.nat);
   printf("Nnres = %d\n", M.nres);

}  /* End of RepeatStruc */

/**********************************************************************
**********************************************************************/
void BuildInfo()

{  printf("\n\tThe builded polysaccharide %s contains:\n", M.titel);
   printf("\n\t%5d residues and\n", M.nres-1);
   printf("\t%5d atoms\n", M.nat);

}  /* End of BuildInfo */

/* End of file */
