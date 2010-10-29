/*
 * This module contains the source code for exporting cartesian
 * coordinates generated by POLYS to a number of different standards
 * Written by Soren Balling Engelsen, INRA-93.
 * Last revision: S.B.E. 1997
 */
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "defs.h"
#include "typs.h"
#include "extern.h"
#include "fileoper.h"
 
/*
 *      Write MONOBANK (binary) coordinate file
 *      =======================================
 *
 * Write cartesian coordinates for the molecule in polys ".x" format.
 * If *fname is NULL pointer the coordinates will be written to "stdout".
 */
void
write_NMB(char *fname)
{
   FILE *fp;
   int   i;

   if (fname != NULL) 
   {  if ((fp = fileopen(fname, "wb")) != NULL)
         printf("MONOBANK file \"%s\" opened for write access\n", fname);
      fprintf(fp, "%s\n", M.id);
      fprintf(fp, "%s\n", M.titel);
      fprintf(fp, "%5d %5d\n", M.nat, M.nres-1);
      fwrite((void *) A, sizeof(A[0]), M.nat, fp);
      fileclose(fname, fp);
      printf("%s succesfully stored in MONOBANK file \"%s\"\n", M.id, fname);
   }
}

/*
 *      Write MONOBANK (ASCII) coordinate file
 *      ======================================
 *
 * Write cartesian coordinates for the molecule in polys ".x" format.
 * If *fname is NULL pointer the coordinates will be written to "stdout".
 */
void
write_MB(char *fname)
{
   FILE *fp;
   int   i;

   if (fname != NULL)
   {  if ((fp = fileopen(fname, "wt")) != NULL)
         printf("MONOBANK file \"%s\" opened for write access\n", fname);
   }
   else
      fp = stdout;
   fprintf(fp, "%s\n", M.id);
   fprintf(fp, "%s\n", M.titel);
   fprintf(fp, "%5d%5d\n", M.nat, M.nres-1);
   for ( i=0; i<M.nat; i++)
      fprintf(fp, "%5d %3d %15.5f %15.5f %15.5f %7.4f %6d %s\n",
                  i, A[i].res, A[i].pos.x, A[i].pos.y, A[i].pos.z,
                     A[i].chg, A[i].ztyp, A[i].lab);
   if (fp != stdout)
   {  fileclose(fname, fp);
      printf("Data for %s succesfully stored in file \"%s\"\n", M.id, fname);
   }
}

/*
 *      Write POLYS format coordinate file
 *      ==================================
 *
 * Write cartesian coordinates for the molecule in polys ".x" format.
 * If *fname is NULL pointer the coordinates will be written to "stdout".
 */
void
write_PF(char *fname, BOOLEAN descarte)
{
   FILE     *fp;
   BOOLEAN  newline=TRUE;
   int      i, j, k;
   Vector3  av;

   if (fname != NULL) 
   {  if ((fp = fileopen(fname, "wt")) != NULL)
         printf("Molecular file \"%s\" opened for write access\n", fname);
   }
   else
      fp = stdout;
   fprintf(fp, "%s\n", M.id);
   fprintf(fp, "%s\n", M.titel);

   if (descarte == TRUE)
   {  CELL.a = 1.0; CELL.b = 1.0; CELL.c = 1.0;
      CELL.alpha = 90.0; CELL.beta = 90.0; CELL.gamma = 90.0;
   }
   else
   {  printf("Storing as fractional coordinates\n"); 
      InvMat3x3(CTF);
      for (i=0; i<M.nat; i++)
      {  av = A[i].pos;
         TransfV(av, CTF, &A[i].pos);
      }
      printf("Transformation Matrix:\n\n");
      matprint(CTF);
   }

   fprintf(fp, "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n",
               CELL.a, CELL.b, CELL.c, CELL.alpha, CELL.beta, CELL.gamma);
   fprintf(fp, "%5d%5d%5d\n", M.nat, M.nbond, nmat);
   for ( i=0; i<M.nat; i++)
      fprintf(fp, "%5d %3d %15.5f %15.5f %15.5f %7.4f %6d %s\n",
                  i, A[i].res, A[i].pos.x, A[i].pos.y, A[i].pos.z, 
                     A[i].chg, A[i].ztyp, A[i].lab);
   for ( i=0; i<M.nbond; i++)
   { newline = FALSE;
     fprintf(fp, "%5d%5d", B[i].from, B[i].to);
     if (((i+1) % 8) == 0)  
     { newline = TRUE;
       fprintf(fp, "\n");
     }
   }
   if (!newline)  
      fprintf(fp, "\n");

   /* if fractional format write symmetry operations */
   if (descarte == FALSE)
     for ( k=0; k<nmat; k++)
       for ( i=0; i<4; i++ )
       { for ( j=0; j<4; j++ )
           fprintf(fp, "%10.6f", MAT[k][i][j]);
         fprintf(fp, "\n");
       }

   if (fp != stdout)
   {  fileclose(fname, fp);
      printf("Data for %s succesfully stored in file \"%s\"\n", M.id, fname);
   }
} /* End of write_PF */

/* Writes CFF coordinate files */
void
write_CFF(char *fname)
{
   register int i;
   int      atyp, btype;
   char     styp[4];
   FILE     *fp;
   
   if (fname != NULL) 
   {  if ((fp = fileopen(fname, "wt")) != NULL)
         printf("Molecular file \"%s\" opened for write access\n", fname);
   }
   else
      fp = stdout;

   fprintf(fp, "%s\n", M.titel);
   fprintf(fp, "F O R M U L A\n");
   fprintf(fp, " 1       0       099 0 1 1 1\n");

   for ( i=0; i<M.nat; i++)
   {  atyp = (A[i].ztyp / 100);
      fprintf(fp, "%20.10f%20.10f%20.10f %3d %-10s\n",
                   A[i].pos.x, A[i].pos.y, A[i].pos.z, atyp, A[i].lab);
   }

   fprintf(fp, "End of >%s<\n", M.id);

   if (fp != stdout)
   {  fileclose(fname, fp);
      printf("Data for %s succesfully stored in file \"%s\"\n", M.id, fname);
   }
}

/*
 *      Writes DAT format coordinate files.
 *      =====================================
 *
 * Write cartesian coordinates for the molecule in CRYSTAL ".dat" format.
 * If *fname is NULL pointer the coordinates will be written to "stdout".
 *
 * 1 LINE:  title
 * 2 LINE:  orthogonal transformations
 * 3 LINE:  lines per atom
 * LINE:  A[i].lab, A[i].pos.x, A[i].pos.y, A[i].pos.z, type, A[i].chg, res
 */
void
write_DAT(char *fname, BOOLEAN descarte)
{
   register int i;
   int      atyp;
   Vector3  av;
   FILE     *fp;
   
   if (fname != NULL) 
   {  if ((fp = fileopen(fname, "wt")) != NULL)
         printf("Molecular file \"%s\" opened for write access\n", fname);
   }
   else
      fp = stdout;

   fprintf(fp, "%s\n", M.titel);

   if (descarte == TRUE)
   {  CELL.a = 1.0; CELL.b = 1.0; CELL.c = 1.0;
      CELL.alpha = 90.0; CELL.beta = 90.0; CELL.gamma = 90.0;
   }
   else
   {  printf("Storing as fractional coordinates\n"); 
      InvMat3x3(CTF);
      for (i=0; i<M.nat; i++)
      {  av = A[i].pos;
         TransfV(av, CTF, &A[i].pos);
      }
      printf("Transformation Matrix:\n\n");
      matprint(CTF);
   }

   fprintf(fp, "%9.4f,%9.5f,%9.5f,%9.5f,%9.5f,%9.5f\n",
               CELL.a, CELL.b, CELL.c, CELL.alpha, CELL.beta, CELL.gamma);
   fprintf(fp, " 1\n");
   for ( i=0; i<M.nat; i++)
   {  atyp = A[i].ztyp;
      switch(atyp)
      /*                                  element            */
      {  case  101: atyp =  3; break;  /* aliphatic hydrogen */
         case  111: atyp = -3; break;  /* hydroxyl hydrogen  */
         case  603: atyp =  1; break;  /* carbon sp2         */
         case  604: atyp =  1; break;  /* carbon sp3         */
         case  614: atyp =  1; break;  /* anomeric carbon    */
         case  801: atyp =  4; break;  /* carbonyl oxygen    */
         case  802: atyp =  2; break;  /* hydroxyl oxygen    */
         case  812: atyp =  4; break;  /* ether oxygen       */
         case  822: atyp =  4; break;  /* carboxyl oxygen    */
         case  703: atyp =  6; break;  /* nitrogen           */
         case 1602: atyp =  5; break;  /* sulphur sp2        */
         case 1604: atyp =  5; break;  /* sulphur sp3        */
         default:   atyp =  9; break;  /* dummy              */
      }
      fprintf(fp, " %-8s,%9.5f,%9.5f,%9.5f,%6d,%9.5f,%5d\n",
                    A[i].lab, A[i].pos.x, A[i].pos.y, A[i].pos.z, atyp, 
                    A[i].chg, A[i].res);
   }
   fprintf(fp, " fin     ,%9.5f,%9.5f,%9.5f,%6d,%9.5f,%5d\n",
              A[i].pos.x, A[i].pos.y, A[i].pos.z, atyp, 
              A[i].chg, A[i].res);
   if (fp != stdout)
   {  fileclose(fname, fp);
      printf("Data for %s succesfully stored in file \"%s\"\n", M.id, fname);
   }
} /* End of write_DAT */

/*
 *      Writes CHARMM format coordinate files.
 *      =======================================
 * Write cartesian coordinates for the molecule in CHARMM ".CRD" format.
 * If *fname is NULL pointer the coordinates will be written to "stdout".
 */
void
write_CHARMM(char *fname)
{
   register int i;
   int      yy, mm, dd, hh, min, ss;
   char     styp[4];
   FILE     *fp;
   
   if (fname != NULL) 
   {  if ((fp = fileopen(fname, "wt")) != NULL)
         printf("Molecular file \"%s\" opened for write access\n", fname);
   }
   else
      fp = stdout;

   fprintf(fp, "*\tTitel:\t\t\t%s\n", M.titel);
   fprintf(fp, "*\tMade by:\t\t%s\n", UserName);
   dtime(&yy,&mm,&dd,&hh,&min,&ss);
   fprintf(fp, "*\tCreation time:\t\t%d:%02d on %d-%d-19%d\n",
                                     hh,min,    dd,mm,  yy);
   fprintf(fp, "*\n");
   fprintf(fp, "%5d\n", M.nat);

   strcpy(styp, "SACC");
   for ( i=0; i<M.nat; i++)
   {  
      fprintf(fp, "%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %d %12.5f\n",
                   i+1, A[i].res, styp, A[i].lab, A[i].pos.x, A[i].pos.y,
                        A[i].pos.z, styp, A[i].res, A[i].chg);
   }

   if (fp != stdout)
   {  fileclose(fname, fp);
      printf("Data for %s succesfully stored in file \"%s\"\n", M.id, fname);
   }
} /* End of write_CHARMM */

/*
 * The MM3 input topology need cycles in its connected lists.
 * This recursive function finds a new ring connection (previous-ai) 
 * from the array cycat[][0] to progress the closure of the ring. 
 * The new ring atom (ai) is checked for covalent bond to the initial 
 * ring atom (first) - if so the ring is closed and the function 
 * exited.
 * BEWARE the search for ring closure is terminated if it has not been
 * found within 16 atoms (which is also the limit of the MM3 connected 
 * lists).
 * BEWARE this function is most likely to fail with systems with double 
 * rings sharing one or more atoms, however the other connected lists should
 * be error free and just some manual adjustment is needed to form the
 *ring closures (bDManPyr, aDAnGalp, cyclo dextrin ....).
 */
int
rec_ring(int cycat[700][2], int connmm3[51][17], int ncon, int count,
         int first, int previous)
{
   int    ia, ncount;
   double dij;

   if (count > 16)
   {  printf("WARNING:\n");
      printf("WARNING: POLYS could not find ring closure\n");
      printf("WARNING: within 16 bonds\n");
      printf("WARNING:\n");
      exit(0);
   }

   for (ia=0; ia<M.nat; ia++)
   {  /* a chain atom and bonded to previous chain atom */
      if ((cycat[ia][0] > 1) AND (actiontype(A[previous],A[ia],&dij) == 1))
      {  connmm3[ncon][count++] = ia;
         cycat[ia][0]--;
         cycat[ia][0]--;
         /* check for ring closure */
         if (actiontype(A[first],A[ia],&dij) == 1) 
         {  connmm3[ncon][count++] = first;  /* special MM3 closure */
            printf("%d membered ring closure\n", count-2);
            return (count);
         }
         ncount = count;
         count = rec_ring(cycat, connmm3, ncon, ncount, first, ia);
         return (count);
      }
   }

   printf("WARNING:\n");
   printf("WARNING: POLYS could not find ring closure\n");
   printf("WARNING:\n");
   return (count);
} /* End of rec_ring */

/*
 * The MM3 input topology need connected lists of chain atoms (more than
 * one bond).
 * This recursive function finds the next connection in the (previous-ia) 
 * from the array cycat[][1] to progress the chain until it is terminated
 * by a terminal atom.
 * The recursive search for terminal atoms is terminated if it has not 
 * been found within 16 atoms (which is also the limit of the MM3 connected
 * lists), the chain will thus be divided into two or more chains.
 */
int
rec_list(int cycat[700][2], int connmm3[51][17], 
         int ncon, int count, int previous)
{
   int    ia, ncount;
   double dij; 

   if (count > 16)
   {  printf("BEWARE:\n");
      printf("BEWARE: MM3 connected list exceeds 16 members,\n");
      printf("BEWARE: it is being divided\n");
      printf("BEWARE:\n");
      return(count);
   }

   for (ia=0; ia<M.nat; ia++)
   {  if ((cycat[ia][1] > 0) AND (ia!=previous))
         if ((actiontype(A[previous],A[ia],&dij) == 1) AND
             ((cycat[previous][0] == 0) OR (cycat[ia][0] == 0)))
      {  /* update connection list */
         cycat[previous][1]--;
         cycat[ia][1]--;
         connmm3[ncon][count++] = ia;
         ncount = count;
         count = rec_list(cycat, connmm3, ncon, ncount, ia);
         return (count);
      }
   }
   return (count);
} /* End of rec_list */

/*
 *      Writes MM3 format coordinate files.
 *      =====================================
 * Write cartesian coordinates for the molecule in MM3 ".MM3" format.
 * If *fname is NULL pointer the coordinates will be written to "stdout".
 *
 * 1* LINE:  (6A10,I1,I4,I2,I1,2I2,I3,F5.0) Name card
 * 2* LINE:  (I1,I4,5X,F10.5,9I5,5X,2I5)    Lists and options
 * 5* LINE:  (3F10.5, I5)                   Coordinates
 *         X coordinate (%10.5lf)
 *         Y coordinate (%10.5lf)
 *         Z coordinate (%10.5lf)
 *         atom type letter (%c)           Imberty
 *         atom type (%5d)
 *         atom label (%s)                 Imberty
 */
void
write_MM3(char *fname)
{
   register int i, j;
   int      i1, i2;
   int      cycle, level, count, ncount;
   int      yy, mm, dd, hh, min, ss;
   int      ncon, ncyc, nattch, attype;
   int      connmm3[51][17], attch[700][3], cycat[700][2];
   double   del, dij;
   FILE     *fp;

   if (fname != NULL) 
   {  if ((fp = fileopen(fname, "wt")) != NULL)
         printf("Molecular file \"%s\" opened for write access\n", fname);
   }
   else
      fp = stdout;

   /*  nattch - bonds to sideatoms            */
   /*  cycat[700][0] == number of ring bonds  */
   /*  cycat[700][1] == number of chain bonds */
   nattch = 1;

   for (i=0; i<M.nbond; i++)
   {  if ((A[B[i].from].ztyp % 10 == 1) OR (A[B[i].to].ztyp % 10 == 1)) 
      {  attch[nattch][1] = B[i].from;
         attch[nattch][2] = B[i].to;
         nattch++;
      }
      else
      {  cycle = FindConn(M.nat, M.nbond, B[i].from, B[i].to, &level);
         /* mark number of ring connections */
         if (cycle == 0)
         {  cycat[B[i].from][0]++;
            cycat[B[i].to][0]++;
         }
         else
         /* or number of chain atom connections */
         {  cycat[B[i].from][1]++;
            cycat[B[i].to][1]++;
         }
      }
   }

   /*  ncon - number of chains of connected atoms */
   ncon = 1; count = 1;

   /* now find the connected lists */
   count = 1;
   for (i1=0; i1<M.nat; i1++)
   {  /* let take a terminal atom */
      if (cycat[i1][1] == 1)
      {  /* update connection list */
         connmm3[ncon][count++] = i1;
         cycat[i1][1]--;
         for (i2=0; i2<M.nat; i2++)
         {  if ((cycat[i2][1] > 0) AND (i2 != i1))
               if ((actiontype(A[i1],A[i2],&dij) == 1) AND
                   ((cycat[i1][0] == 0) OR (cycat[i2][0] == 0)))
            {  /* update connection list */
               connmm3[ncon][count++] = i2;
               cycat[i2][1]--;
               ncount = rec_list(cycat, connmm3, ncon, count, i2);
            }
         }
         connmm3[ncon][0] = ncount;
         ncon++;
         count = 1;
      }
   } 

   /*  recursive search for ring closures */
   for (i1=0; i1<M.nat; i1++)
   {  if (cycat[i1][0] > 1)
      {  connmm3[ncon][count++] = i1;
         cycat[i1][0]--;
         cycat[i1][0]--;
         for (i2=0; i2<M.nat; i2++)
         {  if ((cycat[i2][0] > 0) AND (actiontype(A[i1],A[i2],&dij) == 1))
            {  connmm3[ncon][count++] = i2;
               cycat[i2][0]--;
               cycat[i2][0]--;
               ncount = rec_ring(cycat, connmm3, ncon, count, i1, i2);
            }
         } /* atom 2 */
         connmm3[ncon][0] = ncount;
         ncon++;
         count = 1;
      }
   } /* atom 1 */

   /* if we have superflourous ring connections then convert */

   /* now find the rest of the connected lists */
   count = 1;
   for (i1=0; i1<M.nat; i1++)
   {  /* let take a terminal atom */
      if (cycat[i1][1] == 1)
      {  /* update connection list */
         connmm3[ncon][count++] = i1;
         cycat[i1][1]--;
         for (i2=0; i2<M.nat; i2++)
         {  if ((cycat[i2][1] > 0) AND (i2 != i1))
               if ((actiontype(A[i1],A[i2],&dij) == 1) AND
                   ((cycat[i1][0] == 0) OR (cycat[i2][0] == 0)))
            {  /* update connection list */
               connmm3[ncon][count++] = i2;
               cycat[i2][1]--;
               ncount = rec_list(cycat, connmm3, ncon, count, i2);
            }
         }
         connmm3[ncon][0] = ncount;
         ncon++;
         count = 1;
      }
   } 

/* Modification by M. Kouwijzer: energy will be minimized ( so del
   has been changed from 0.0 to 0.00003) instead of only the
   calculation of the initial energy */

   del = 0.00003;
/* temporary change: only an initial energy is calculated */
   fprintf(fp, "%-60s0%4d 4  0 1  1999.0\n", M.titel, M.nat);
/*   fprintf(fp, "%-60s0%4d 4  0 0  1999.0\n", M.titel, M.nat);
   end temporary change */
   fprintf(fp, "0%4d     %10.5lf     %5d", ncon-1, del, nattch-1);
   fprintf(fp,"    0         0    0    0    0    0         1    0\n");

   for (i=1; i<ncon; i++)
   {  for (j=1; j<connmm3[i][0]; j++)
         fprintf(fp,"%5d", connmm3[i][j]+1);
      fprintf(fp, "\n");
   }

   for ( j=1; j<nattch; j++)
   {  fprintf(fp,"%5d%5d", attch[j][1]+1, attch[j][2]+1);
      if ((j % 8) == 0)
         fprintf(fp, "\n");
   }
   if (((nattch-1) % 8) != 0)
      fprintf(fp, "\n");
 
   /* convert atomtypes on the fly */
   attype = 0;
   for ( i=0; i<M.nat; i++)
   {  
      switch(A[i].ztyp)
      {  case  101: attype=5;  break; /* except on N, O, S */
         case  111: attype=21; break; /* -OH */
         case  121: attype=28; break; /* amide */
         case  131: attype=23; break; /* amine */
         case  141: attype=44; break; /* S-H Thiol */
         case  602: attype=4;  break; /* alkyne */
         case  603: attype=3;  break; /* carbonyl */
         case  613: attype=2;  break; /* alkene */
         case  623: attype=50; break; /* aromatic - benzene */
         case  604: attype=1;  break; /* aliphatic */
         case  614: attype=1;  break; /* aliphatic - anomeric */
         case  701: attype=10; break; /* azo nitrogen */
         case  703: attype=9;  break; /* sp2 nitrogen */
         case  713: attype=39; break; /* amine unprotonated */
         case  723: attype=8;  break; /* sp3 nitrogen */
         case  704: attype=39; break; /* amine */
         case  801: attype=7;  break; /* carbonyl */
         case  811: attype=47; break; /* carboxylate ion */
         case  802: attype=6;  break; /* hydroxyl */
         case  812: attype=6;  break; /* ether */
         case  822: attype=75; break; /* carboxyl */
         case 1604: attype=18; break; /* sulphur */
         case 1701: attype=12; break; /* chloride */
           default: attype=0;  
                    printf("WARNING:\n");
                    printf("WARNING: Undetermined atomtype conversion\n");
                    printf("WARNING:\n");
                    break;
      }
      fprintf(fp, "%10.5lf%10.5lf%10.5lf %4d%s\n",
              A[i].pos.x, A[i].pos.y, A[i].pos.z, attype, A[i].lab);
   }

   if (fp != stdout)
   {  fileclose(fname, fp);
      printf("Data for %s succesfully stored in file \"%s\"\n", M.id, fname);
   }
} /* End of write_MM3 */

/*
 *      Writes SYBYL/PIM format coordinate files.
 *      =========================================
 * Write cartesian coordinates for the molecule in SYBYL ".mol2" format.
 * If *fname is NULL pointer the coordinates will be written to "stdout".
 *
 * --------------------------------------------------
 * Elements   SYBYL/PIM        dat      polys
 * --------------------------------------------------
 * Carbon    sp3 =  C.3/C.C     1         604
 * Carbon    sp3 =  C.3/C.A     1         614 (a-anomeric)
 * Carbon    sp3 =  C.3/C/B     1         614 (b-anomeric)
 * Carbon    sp2 =  C.2         1         603
 * Carbon    sp  =  C.1         1         602
 * Carbon    ar  =  C.ar        1         623
 * --------------------------------------------------
 * Oxygen    sp3 =  O.3         2         802 (hydroxyl)
 * Oxygen    sp3 =  O.3         2         822 (carboxyl-O)
 * Oxygen    sp3 =  O.3/O.A     4         812 (a-anomeric)
 * Oxygen    sp3 =  O.3/O.B     4         812 (b-anomeric)
 * Oxygen    sp3 =  O.3/O.R     4         812 (sugar ring-O)
 * Oxygen    sp2 =  O.2         4         801
 * --------------------------------------------------
 * Hydrogen      =  H/H.C       3         101
 * Hydrogen      =  H          -3         111
 * Hydrogen      =  H           3         121
 * --------------------------------------------------
 * Nitrogen  sp3 =  N.3         5         703
 * Nitrogen  sp2 =  N.2         5         713
 * Nitrogen  sp  =  N.1         5         702
 * Nitrogen  ar  =  N.ar        5         723
 * Nitrogen  amid=  N.am        5         703
 * Nitrogen  plan=  N.pl3       5         723
 * Nitrogen  pos =  N.4         5         704
 * --------------------------------------------------
 * Sulfur    sp3 =  S.3         6        1602
 * Sulfur    sp2 =  S.2         6        1601
 * Sulfur    sp3 =  S.O        18        1604
 * Sulfur    sp3 =  S.O2        6        1603
 * --------------------------------------------------
 * Phosphour sp3 =  P.3         ?        1504
 * --------------------------------------------------
 * Dummy         =  Du          ?           9
 *--------------------------------------------------
 */
void
write_SYBYL(char *fname)
{
   register int i;
   int      atyp, btype;
   int      yy, mm, dd, hh, min, ss;
   char     styp[4];
   FILE     *fp;
   
   if (fname != NULL) 
   {  if ((fp = fileopen(fname, "wt")) != NULL)
         printf("Molecular file \"%s\" opened for write access\n", fname);
   }
   else
      fp = stdout;

   fprintf(fp, "#\tName:\t\t\t%s\n", M.titel);
   fprintf(fp, "#\tCreating user name:\t%s\n", UserName);
   dtime(&yy,&mm,&dd,&hh,&min,&ss);
   fprintf(fp, "#\tCreation time:\t\t%d:%02d on %d-%d-%d\n\n",
                                     hh,min,    dd,mm,  yy);
   fprintf(fp, "#\tModifying user name:\t%s\n", UserName);
   fprintf(fp, "#\tModification time:\t%d:%02d on %d-%d-%d\n\n",
                                         hh,min,    dd,mm,  yy);

   fprintf(fp, "@<TRIPOS>MOLECULE\n");
   fprintf(fp, "%s\n", M.id);
   fprintf(fp, "%5d %5d %5d %5d %5d\n", M.nat, M.nbond, 1, 0, 0);
   fprintf(fp, "SMALL\n");
   fprintf(fp, "USER_CHARGES\n\n\n");

   fprintf(fp, "@<TRIPOS>ATOM\n");
   for ( i=0; i<M.nat; i++)
   {  atyp = A[i].ztyp;
      switch(atyp)
      {  case  101: strcpy(styp, "H   "); break;  /* hydrogen     */
         case  111: strcpy(styp, "H   "); break;  /* hydrogen     */
         case  121: strcpy(styp, "H   "); break;  /* hydrogen     */
         case  109: strcpy(styp, "H   "); break;  /* hydrogen     */
         case  603: strcpy(styp, "C.2 "); break;  /* sp2 carbon   */
         case  604: strcpy(styp, "C.3 "); break;  /* sp3 carbon   */
         case  614: strcpy(styp, "C.3 "); break;  /* a-anomeric C */
         case  703: strcpy(styp, "N.am"); break;  /* sp3 nitrogen */
         case  801: strcpy(styp, "O.2 "); break;  /* sp2 oxygen   */
         case  811: strcpy(styp, "O.2 "); break;  /* sp2 oxygen   */
         case  802: strcpy(styp, "O.3 "); break;  /* sp3 oxygen   */
         case  812: strcpy(styp, "O.3 "); break;  /* sugar ring-O */
         case  822: strcpy(styp, "O.3 "); break;  /* carboxyl-O */
         case 1602: strcpy(styp, "S.2 "); break;  /* sp2 sulphur  */
         case 1604: strcpy(styp, "S.O "); break;  /* sp2 sulphur  */
         default:   strcpy(styp, "DU  "); break;  /* dummy        */
      }
      fprintf(fp, "%7d %-8s %9.4f %9.4f %9.4f %-4s %6d %-8s %8.4f\n",
                   i+1, A[i].lab, A[i].pos.x, A[i].pos.y, A[i].pos.z, styp, 
                   (A[i].res+1), "<1>", A[i].chg);
   }

   fprintf(fp, "@<TRIPOS>BOND\n");
   /* sybyl atom numbering offset by one */
   btype = 1;          /* bond order  */
   for ( i=0; i<M.nbond; i++)
   {  if ((A[B[i].from].ztyp == 801) OR (A[B[i].to].ztyp == 801)
         OR (A[B[i].to].ztyp == 811) OR (A[B[i].to].ztyp == 811))
         btype = 2;
      fprintf(fp, "%6d %4d %4d %d\n", i+1, B[i].from+1, B[i].to+1, btype);
      btype = 1;
   }

   fprintf(fp, "@<TRIPOS>SUBSTRUCTURE\n");
   fprintf(fp, "     1 NO_NAME     1 PERM              0 ****  ****    0 ROOT\n");

   if (fp != stdout)
   {  fileclose(fname, fp);
      printf("Data for %s succesfully stored in file \"%s\"\n", M.id, fname);
   }
} /* End of write_SYBYL */

/* */
void
write_PIM(char *fname)
{
   register int i;
   int      atyp, btype;
   int      yy, mm, dd, hh, min, ss;
   char     styp[4];
   FILE     *fp;

   if (fname != NULL)
   {  if ((fp = fileopen(fname, "wt")) != NULL)
         printf("Molecular file \"%s\" opened for write access\n", fname);
   }
   else
      fp = stdout;

   fprintf(fp, "#\tName:\t\t\t%s\n", M.titel);
   fprintf(fp, "#\tCreating user name:\t%s\n", UserName);
   dtime(&yy,&mm,&dd,&hh,&min,&ss);
   fprintf(fp, "#\tCreation time:\t\t%d:%02d on %d-%d-%d\n\n",
                                     hh,min,    dd,mm,  yy);
   fprintf(fp, "#\tModifying user name:\t%s\n", UserName);
   fprintf(fp, "#\tModification time:\t%d:%02d on %d-%d-%d\n\n",
                                         hh,min,    dd,mm,  yy);

   fprintf(fp, "@<TRIPOS>MOLECULE\n");
   fprintf(fp, "%s\n", M.id);
   fprintf(fp, "%5d %5d %5d %5d %5d\n", M.nat, M.nbond, 1, 0, 0);
   fprintf(fp, "SMALL\n");
   fprintf(fp, "USER_CHARGES\n\n\n");

   fprintf(fp, "@<TRIPOS>ATOM\n");
   for ( i=0; i<M.nat; i++)
   {  atyp = A[i].ztyp;
      switch(atyp)
      /*
       * ATOMTYPES converted to PIM types (NOT standard SYBYL)
       * Imberty, Hardman, Carver and Perez.
       * Glycobiology, 1 (6), pp. 631-42, 1991.
       * BEWARE. Anomeric carbons are always converted to C.A
       * BEWARE. Ether oxygens are always converted to O.R 
       */
      {  case  101: strcpy(styp, "H.C "); break;  /* hydrogen     */
         case  111: strcpy(styp, "H   "); break;  /* hydrogen     */
         case  121: strcpy(styp, "H   "); break;  /* hydrogen     */
         case  109: strcpy(styp, "H   "); break;  /* hydrogen     */
         case  603: strcpy(styp, "C.2 "); break;  /* sp2 carbon   */
         case  604: strcpy(styp, "C.C "); break;  /* sp3 carbon   */
         case  614: strcpy(styp, "C.A "); break;  /* a-anomeric C */
         case  703: strcpy(styp, "N.am"); break;  /* sp3 nitrogen */
         case  801: strcpy(styp, "O.2 "); break;  /* sp2 oxygen   */
         case  811: strcpy(styp, "O.2 "); break;  /* sp2 oxygen   */
         case  802: strcpy(styp, "O.3 "); break;  /* sp3 oxygen   */
         case  812: strcpy(styp, "O.R "); break;  /* sugar ring-O */
         case  822: strcpy(styp, "O.3 "); break;  /* carboxyl-O */
         case 1602: strcpy(styp, "S.2 "); break;  /* sp2 sulphur  */
         case 1604: strcpy(styp, "S.O "); break;  /* sp2 sulphur  */
         default:   strcpy(styp, "DU  "); break;  /* dummy        */
      }
      fprintf(fp, "%7d %-8s %9.4f %9.4f %9.4f %-4s %6d %-8s %8.4f\n",
                   i+1, A[i].lab, A[i].pos.x, A[i].pos.y, A[i].pos.z, styp,
                   (A[i].res+1), "<1>", A[i].chg);
   }

   fprintf(fp, "@<TRIPOS>BOND\n");
   /* sybyl atom numbering offset by one */
   btype = 1;          /* bond order  */
   for ( i=0; i<M.nbond; i++)
   {  if ((A[B[i].from].ztyp == 801) OR (A[B[i].to].ztyp == 801)
         OR (A[B[i].to].ztyp == 811) OR (A[B[i].to].ztyp == 811))
         btype = 2;
      fprintf(fp, "%6d %4d %4d %d\n", i+1, B[i].from+1, B[i].to+1, btype);
      btype = 1;
   }

   fprintf(fp, "@<TRIPOS>SUBSTRUCTURE\n");
   fprintf(fp, "     1 NO_NAME     1 PERM              0 ****  ****    0 ROOT\n");

   if (fp != stdout)
   {  fileclose(fname, fp);
      printf("Data for %s succesfully stored in file \"%s\"\n", M.id, fname);
   }
} /* End of write_PIM */

/* Writes TINKER XYZ format coordinate files */
void
write_XYZ(char *fname)
{
   register int i,j;
   int      attype;
   int      yy, mm, dd, hh, min, ss;
   FILE     *fp;
   
   if (fname != NULL) 
   {  if ((fp = fileopen(fname, "wt")) != NULL)
         printf("Molecular file \"%s\" opened for write access\n", fname);
   }
   else
      fp = stdout;

   fprintf(fp, "%6d %s\n", M.nat, M.titel);
   for ( i=0; i<M.nat; i++)
   {  switch(A[i].ztyp)
      {  case  101: attype=5;  break; /* except on N, O, S */
         case  111: attype=21; break; /* -OH */
         case  121: attype=28; break; /* amide */
         case  131: attype=23; break; /* amine */
         case  141: attype=44; break; /* S-H Thiol */
         case  602: attype=4;  break; /* alkyne */
         case  603: attype=3;  break; /* carbonyl */
         case  613: attype=2;  break; /* alkene */
         case  623: attype=50; break; /* aromatic - benzene */
         case  604: attype=1;  break; /* aliphatic */
         case  614: attype=1;  break; /* aliphatic - anomeric */
         case  701: attype=10; break; /* azo nitrogen */
         case  703: attype=9;  break; /* sp2 nitrogen */
         case  713: attype=39; break; /* amine unprotonated */
         case  723: attype=8;  break; /* sp3 nitrogen */
         case  704: attype=39; break; /* amine */
         case  801: attype=7;  break; /* carbonyl */
         case  811: attype=47; break; /* carboxylate ion */
         case  802: attype=6;  break; /* hydroxyl */
         case  812: attype=6;  break; /* ether */
         case  822: attype=75; break; /* carboxyl */
         case 1604: attype=18; break; /* sulphur */
         case 1701: attype=12; break; /* chloride */
           default: attype=0;  
                    printf("WARNING:\n");
                    printf("WARNING: Undetermined atomtype conversion\n");
                    printf("WARNING:\n");
                    break;
      }
      /* Peter find en loesning paa udskrivning af %s-3 label */
      fprintf(fp, "%6d  %c%c%c%12.6f%12.6f%12.6f%6d",
              i+1, A[i].lab[0], A[i].lab[1], A[i].lab[2],
                   A[i].pos.x, A[i].pos.y, A[i].pos.z, attype);

      for (j=0; j<M.nbond; j++)
      {  if (B[j].from == i)
         {
            fprintf(fp, "%6d", B[j].to);
         }
         if (B[j].to == i) 
         {
            fprintf(fp, "%6d", B[j].from);
         }
      }
      fprintf(fp, "\n");

   }

   if (fp != stdout)
   {  fileclose(fname, fp);
      printf("Data for %s succesfully stored in file \"%s\"\n", M.id, fname);
   }
} /* End of write_XYZ */

/*
 *       Writes SPF format coordinate files.
 *      =======================================
 * Write cartesian coordinates for the molecule in ".spf" format, to 
 * be used with the programs PLUTON and PLATON by Anthony L. Spek 
 * (A.L. Spek, Acta Cryst. A46 (1990) C34).
 * If *fname is NULL pointer the coordinates will be written to "stdout".
 */
void
write_SPF(char *fname)
{
   register int i;
   int      n, nh=0, nc=0, nn=0, no=0;
   FILE     *fp;
   
   if (fname != NULL) 
   {  if ((fp = fileopen(fname, "wt")) != NULL)
         printf("Molecular file \"%s\" opened for write access\n", fname);
   }
   else
      fp = stdout;

   fprintf(fp, "TITLE %s\n", M.titel);
   fprintf(fp, "CELL 999.0 999.0 999.0 90.0 90.0 90.0\n");
   fprintf(fp, "SPGR P1\n");
   fprintf(fp, "ANGSTROM\n");

   for ( i=0; i<M.nat; i++)
   {  n = (int) (A[i].ztyp / 100);
      switch(n)
      {  case  1: /* hydrogen */
         nh++;
         fprintf(fp, "ATOM H%-3d %10.5f%10.5f%10.5f\n", nh, 
                 A[i].pos.x, A[i].pos.y, A[i].pos.z);
         break;  
         case  6: /* carbon */
         nc++;
         fprintf(fp, "ATOM C%-3d %10.5f%10.5f%10.5f\n", nc, 
                 A[i].pos.x, A[i].pos.y, A[i].pos.z);
         break;  
         case  7: /* nitrogen */
         nn++;
         fprintf(fp, "ATOM N%-3d %10.5f%10.5f%10.5f\n", nn, 
                 A[i].pos.x, A[i].pos.y, A[i].pos.z);
         break;  
         case  8: /* oxygen */
         no++;
         fprintf(fp, "ATOM O%-3d %10.5f%10.5f%10.5f\n", no, 
                 A[i].pos.x, A[i].pos.y, A[i].pos.z);
         break;  
         default: /* not C, O, N, or H */
         fprintf(fp, "ATOM %-4s %10.5f%10.5f%10.5f\n", A[i].lab, 
                 A[i].pos.x, A[i].pos.y, A[i].pos.z);
         break;  
       }
    }

/*   for ( i=0; i<M.nat; i++)
   {  
      fprintf(fp, "ATOM %-4s %10.5f%10.5f%10.5f\n", A[i].lab, 
            A[i].pos.x, A[i].pos.y, A[i].pos.z);
   } */
   fprintf(fp, "ENDS\n");
   fprintf(fp, "STOP\n");

   if (fp != stdout)
   {  fileclose(fname, fp);
      printf("Data for %s succesfully stored in file \"%s\"\n", M.id, fname);
   }
} /* End of write_SPF */

/*
 *        Writes PDB format coordinate files.
 *      =========================================
 * Write cartesian coordinates for the molecule in PDB format.
 * If *fname is NULL pointer the coordinates will be written to "stdout".
 *
 * MTL 27.01.98 quick and dirty PDB output (good enough for RASMOL)
 */
void
write_PDB(char *fname)
{
  int	i;
  int	atyp, btype;
  int	yy, mm, dd, hh, min, ss;
  char	styp[128], buf1[128];
  FILE	*fp=NULL;

#ifdef toto
   if (fname != NULL)
   {  if ((fp = fileopen(fname, "wt")) != NULL)
         printf("Molecular file \"%s\" opened for write access\n", fname);
   }
   else
      fp = stdout;
#endif

  if(fname!=NULL){
    if((fp=fileopen(fname,"wt"))!=NULL){
      printf("Molecular file \"%s\" opened for write access\n", fname);
    }
  }
  else{
    fp=stdout;
  }
printf("hello mum%d\n",M.nat);
  fprintf(fp,"HEADER    %s\n",M.titel);
  dtime(&yy,&mm,&dd,&hh,&min,&ss);
  fprintf(fp,"AUTHOR    Generated by POLYS %d:%02d on %d-%d-19%d\n",hh,min,dd,mm,yy);

  for(i=0;i<M.nat;i++){
    atyp=A[i].ztyp;
    switch(atyp/100){
      case  1: strcpy(styp, "H"); break;  /* hydrogen     */
      case  6: strcpy(styp, "C"); break;  /* carbon   */
      case  7: strcpy(styp, "N"); break;  /* nitrogen */
      case  8: strcpy(styp, "O"); break;  /* oxygen   */
      case 15: strcpy(styp, "P"); break;  /* phosphorus */
      case 16: strcpy(styp, "S"); break;  /* sulphur  */
      case 17: strcpy(styp, "Cl"); break;  /* sulphur  */
      default:   strcpy(styp, "C"); break;  /* dummy        */
    }

    sprintf(buf1,"%d",i+1);
    strcat(styp,buf1);
    fprintf(fp,"ATOM  %5d %-4.4s UNK  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
        i+1,styp/*A[i].lab*/,A[i].res+1, A[i].pos.x,A[i].pos.y,A[i].pos.z);
  }

/*     fprintf(fp, "%7d %-8s %9.4f %9.4f %9.4f %-4s %6d %-8s %8.4f\n",
                  i+1, A[i].lab, A[i].pos.x, A[i].pos.y, A[i].pos.z, styp,
                  (A[i].res+1), "<1>", A[i].chg);
*/

#ifdef not_for_now
   fprintf(fp, "@<TRIPOS>BOND\n");
  /* sybyl atom numbering offset by one */
  btype = 1;          /* bond order  */
  for(i=0;i<M.nbond;i++){
    if ((A[B[i].from].ztyp == 801) OR (A[B[i].to].ztyp == 801)
         OR (A[B[i].to].ztyp == 811) OR (A[B[i].to].ztyp == 811))
         btype = 2;
      fprintf(fp, "%6d %4d %4d %d\n", i+1, B[i].from+1, B[i].to+1, btype);
      btype = 1;
  }
#endif

  if(fp!=stdout){
    fileclose(fname, fp);
    printf("Data for %s succesfully stored in file \"%s\"\n", M.id, fname);
  }
} /* end of write_PDB() */

/* Writes POLYS plot files */
void
write_PLOT(char *fname)
{
   register int i;
   int      atyp, btype;
   char     styp[4];
   FILE     *fp;
   
   if (fname != NULL) 
   {  if ((fp = fileopen(fname, "wt")) != NULL)
         printf("Molecular file \"%s\" opened for write access\n", fname);
   }
   else
      fp = stdout;

   fprintf(fp, "%5d%5d\n", M.nat, M.nbond);

   for ( i=0; i<M.nat; i++)
   {  atyp = (A[i].ztyp / 100);
      fprintf(fp, "%10.5f%10.5f%10.5f%10d  %-10s\n",
                   A[i].pos.x, A[i].pos.y, A[i].pos.z, atyp, A[i].lab);
   }

   /* sybyl atom numbering offset by one */
   btype = 1;          /* bond order  */
   for ( i=0; i<M.nbond; i++)
   {  if ((A[B[i].from].ztyp == 801) OR (A[B[i].to].ztyp == 801))
         btype = 2;
      fprintf(fp, "%5d%5d%5d\n", B[i].from+1, B[i].to+1, btype);
      btype = 1;
   }

   if (fp != stdout)
   {  fileclose(fname, fp);
      printf("Data for %s succesfully stored in file \"%s\"\n", M.id, fname);
   }
} /* End of write_PLOT */


/* Guides the output rutines */
void
WritCoord(int format, char *fname, BOOLEAN descarte)
{
   /* first remove the quotes from the file name */
   if (fname != NULL)
      DelQuotes(fname);

   printf("Writing Cartesian Coordinates using ");
   switch (format)
   {  case 0     : printf("AMBER Format\n");
                   printf("THIS OPTION IS NOT IMPLEMENTED YET\n");
                   break;
      case 1     : printf("CHARMM Format\n");
                   write_CHARMM(fname); 
                   break;
      case 2     : printf("CFF Format\n");
                   write_CFF(fname);
                   break;
      case 3     : printf("CSD Format\n");
                   printf("CANNOT export to CSD format\n");
                   break;
      case 4     : printf("DAT Format\n");
                   write_DAT(fname, descarte); 
                   break;
      case 5     : printf("MM3 Format\n");
                   write_MM3(fname); 
                   break;
      case 6     : printf("MONOBANK Format\n");
                   write_MB(fname); 
                   break;
      case 7     : printf("PDB Format\n");
                   write_PDB(fname);
                   break;
      case 8     : printf("PF Format\n");
                   write_PF(fname, descarte); 
                   break;
      case 9     : printf("PLOT Format\n");
                   write_PLOT(fname); 
                   break;
      case 10    : printf("SYBYL Format\n");
                   write_SYBYL(fname); 
                   break;
      case 11    : printf("SPF Format\n");
                   write_SPF(fname); 
                   break;
      case 12    : printf("SYBYL PIM Format\n");
                   write_PIM(fname); 
                   break;
      case 13    : printf("TINKER Format\n");
                   write_XYZ(fname); 
   }
}
