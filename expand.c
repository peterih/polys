/**********************************************************************
POLYS MODULE <expand.c>

  This module contains the source code for maintaining and updating 
the residue list Res[] during the PRIMARY command and the code for
3D generation and linking of the primary sequence during the BUILD
command.

  Written by Soren Balling Engelsen, INRA 1994.
  Last revision: S.B.E. 1997
**********************************************************************/
#include <math.h>
#include <string.h>
#include "polys.h"
#include "extern.h"
 


/**********************************************************************
This routine will read a macro structure from a file and attach it to
the POLYSaccharide which has been buildt.
**********************************************************************/
void addmacro(int polyres, int con_b, int con_f, int macrores, char *mf,
              double phi, double psi, double omega)

{  int       i, j, mres, redmax;
   int       sideno[10], next[10];

   printf("Adding a macro structure to the generated polysaccharide..\n");
   printf("Polysaccharide anchoring point:  residue %d position %d\n", 
            polyres, con_b);
printf("Macrostructure attachment site:  residue %d position %d from file %s\n",
          macrores, con_f, mf);

   /* 1. Read the MACRO coordinates in MONOBANK format */
   mres=read_MB(mf, FALSE);

   /* 2. Form the glycosidic linkage */
   GlycLink(M.nres-mres, M.nres-1-mres+macrores, M.nres-1, 
            1,           polyres,                M.nres-mres-1,
            redmax, con_f, con_b, phi, psi, omega); 

   /* 3. Hydrolysis */
   Hydrolysis(da);
 
   /* 4. Reset topology information */
   M.nbond=0;

} /* End of tertiary */

/**********************************************************************
This routine BUILDs the 3D structure of the polysaccharide from the
monosaccharide conectivity table Res[]

Hydrolysis has to be performed seperately by the routine Hydrolysis in
the module <glyclink.c>
**********************************************************************/
void tertiary(int resmax)

{  int       i, j, con_to, oi, dummy;
   int       sideno[10], next[10];

   /* read in the coordinates from MONOBANK */
   for (i=1; i<=resmax; i++)
     dummy=read_MB(Res[i].resid, TRUE);
  
   for (i=0; i<10; i++)
   { sideno[i]=0; 
     next[i]=0;
   }

   for (i=1; i<=resmax; i++)
   { con_to=Res[i].conres;
     if ((Res[i].from != 0))
     { 
       if (debug) 
         printf("\nTERTIARY: --> LINK RESIDUE #%d [S#%d]\n",
                i,sideno[Res[i].status]);

       /* connect after inbetween structure */
       if (next[0] == i)
       {  GlycLink(1, oi, i-1, i, i, i, resmax,
             Res[oi].from, Res[oi].to, Res[oi].phi, Res[oi].psi, Res[oi].omega);
          if (debug) printf("TERTIARY: Connect after inbetween\n");

          for (j=Res[oi].status; j>=0; j--)
            sideno[j]++;
       }

       /* chain prolongation */
       if ((Res[i].status==Res[con_to].status) 
            AND (con_to > i) AND (con_to<=resmax))
       {  if (con_to == i+1)
          {  GlycLink(i-sideno[Res[i].status], i, i, 
                con_to, con_to, con_to, resmax, 
                Res[i].from, Res[i].to, Res[i].phi, Res[i].psi, Res[i].omega);
             for (j=Res[i].status; j>=0; j--)
               sideno[j]++;
          }
          else
          /* wait until the inbetween structure has been built */
          {  next[0] = con_to;
             if (debug) printf("TERTIARY: Inbetween\n");
             oi = i;
          }
       }

       /* side chain attachment */
       if ((Res[i].status > Res[con_to].status) OR (con_to < i))
       {  GlycLink(i-sideno[Res[i].status], i, i, 
                    1, con_to, i-sideno[Res[i].status]-1, resmax,
               Res[i].from, Res[i].to, Res[i].phi, Res[i].psi, Res[i].omega);
          sideno[Res[i].status]=0;
          for (j=Res[i].status-1; j>=0; j--)
            sideno[j]++;
       }
     }
   }

} /* End of tertiary */

/**********************************************************************
This routine modifies connectivity values of individual members 
of the monosaccharide connectivity table
**********************************************************************/
void modify_val(int res, int back, int forward, double theta, 
                double phi, double psi, double omega)

{  int      i;

   Res[res].theta=theta;
   Res[res].phi=phi;
   Res[res].psi=psi;
   Res[res].omega=omega;
   Res[res].from=back;
   Res[res].to=forward;

} /* End of modify_val */

/**********************************************************************
This routine add a new members to the monosaccharide connectivity table

     status:     0   main chain
                 1   side chain
                 2   connecting residue on side chain
**********************************************************************/
void expand(int newres, char *fname, int back, int forward, 
            double theta, double phi, double psi, double omega, 
            int conres, int status)

{  int      i;

   strcpy(Res[newres].resid,fname);
   Res[newres].ringtype=6;            /* pyranose */
   Res[newres].status=status;
   Res[newres].from=back;
   Res[newres].to=forward;
   Res[newres].conres=conres;         /* forward connecting residue */
   Res[newres].theta=theta;
   Res[newres].phi=phi;
   Res[newres].psi=psi;
   Res[newres].omega=omega;

   /* detect if a furanose residue */
   if (strchr(fname, 'f') != NULL)   Res[newres].ringtype=5;
   /* heptose exception for Kdo */
   if (strstr(fname, "Kdo") != NULL) Res[newres].ringtype=7;
   /* heptose exception for Dha */
   if (strstr(fname, "Dha") != NULL) Res[newres].ringtype=7;

   /* save glyclink information in index zero */
   modify_val(0,back,forward,theta,phi,psi,omega);

} /* End of expand */

/**********************************************************************
This routine modifies connectivity parameters of individual members 
of the monosaccharide connectivity table
**********************************************************************/
void modify_con(int res, int conres, int status)

{  int      i;

   if (debug) 
     printf("MODIFY: res #%d from conres [%d] to conres [%d]\n",
             res,Res[res].conres,conres);
   Res[res].conres=conres;
   Res[res].status=status;

} /* End of modify_res */

/**********************************************************************
This routine copies one member of the monosaccharide connectivity table 
Res[source] to another member Res[destin]
**********************************************************************/
void copy_res(int source, int destin, int offset)

{  int      i;

   strcpy(Res[destin].resid,Res[source].resid);
   Res[destin].ringtype= Res[source].ringtype;
   Res[destin].status  = Res[source].status;
   Res[destin].from    = Res[source].from;
   Res[destin].to      = Res[source].to;
   Res[destin].conres  = Res[source].conres+offset;
   Res[destin].theta   = Res[source].theta;
   Res[destin].phi     = Res[source].phi;
   Res[destin].psi     = Res[source].psi;
   Res[destin].omega   = Res[source].omega;

} /* End of copy_res */

/**********************************************************************
This routine prints the monosaccharide connectivity table Res[] on the
standard output
**********************************************************************/
void res_list(int maxres)

{  int      i;

   printf("\nMONOSACCHARIDE CONNECTIVITY TABLE\n\n");
   for (i=1; i<=maxres; i++)
   {  printf("R[%3d]",i);
      if (Res[i].ringtype==7) printf(" H");
      if (Res[i].ringtype==6) printf(" P");
      if (Res[i].ringtype==5) printf(" F");
      printf("%d %-20s (%d->%d) R[%3d] : %7.1lf%7.1lf%7.1lf\n",
              Res[i].status,Res[i].resid,Res[i].from,Res[i].to,Res[i].conres,
              Res[i].phi,Res[i].psi,Res[i].omega);
   }

} /* End of res_list */

/* End of file */
