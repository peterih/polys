/**********************************************************************
  "helix.c"

  This file contains the source code for helix calculations

  Written by Soren Balling Engelsen, INRA-95.
**********************************************************************/
#include <stdio.h>
#include <math.h>
#include "polys.h"
#include "extern.h"


/**********************************************************************
   Calculate the normalized difference vector 
**********************************************************************/
void uvec(Vector3 a, Vector3 b, Vector3 *c)

{  double  norm;

   c->x = a.x - b.x;
   c->y = a.y - b.y;
   c->z = a.z - b.z;
   norm = normalize(&(*c));

} /* End of uvec */

/**********************************************************************
	Calculate the angle: b-a-c
**********************************************************************/
double delta(Vector3 a, Vector3 b, Vector3 c)

{  double   csth, snth;
   Vector3  e, f;

   uvec(a,b,&e);
   uvec(a,c,&f);
   csth=dotprod(e,f);
   snth=sqrt(1.0 - csth*csth);
   if (csth < 0.0)
     return ( atan(snth/csth)*RADtoDEG + 180.0 );
   if (csth == 0.0)
     return ( 90.0 );
   if (csth > 0.0)
     return ( atan(snth/csth)*RADtoDEG );

} /* End of delta */

/**********************************************************************
**********************************************************************/
void mxrot(Vector3 v1, double phi, Matrix Q)

{  double  a, sf, b, c, d, aa, bb, cc, dd, ab, ac, ad, bc, cd, bd;

   a  = cos(phi);
   sf = sin(phi);
   b  = v1.x * sf;
   c  = v1.y * sf;
   d  = v1.z * sf;
   aa = a * a;
   bb = b * b;
   cc = c * c;
   dd = d * d;
   ab = a * b;
   ac = a * c;
   ad = a * d;
   bc = b * c;
   cd = c * d;
   bd = b * d;
   Q[0][0] = aa + bb - cc - dd;
   Q[1][1] = aa - bb + cc - dd;
   Q[2][2] = aa - bb - cc + dd;
   Q[0][1] = 2.0 * (bc - ad);
   Q[1][0] = 2.0 * (bc + ad);
   Q[2][0] = 2.0 * (bd - ac);
   Q[0][2] = 2.0 * (bd + ac);
   Q[1][2] = 2.0 * (cd - ab);
   Q[2][1] = 2.0 * (cd + ab);

} /* End of mxrot */

/**********************************************************************
**********************************************************************/
void convert(Vector3 ai[], int natom, int norg, int nx, int nxy, Vector3 ac[])

{  int       i;
   double    norm, a, b, c;
   Vector3   a1, b1, a2, a3;

   for (i=0; i<natom; i++)
   {  ac[i].x = ai[i].x - ai[norg].x;
      ac[i].y = ai[i].y - ai[norg].y;
      ac[i].z = ai[i].z - ai[norg].z;
   }

   a1   = ac[nx];
   norm = normalize(&a1);
   b1   = ac[nxy];
   norm = normalize(&b1);
   crossprod(a1, b1, &a3);
   norm = normalize(&a3);
   crossprod(a3, a1, &a2);
   norm = normalize(&a2);

   for (i=0; i<natom; i++)
   {  a = dotprod(a1, ac[i]);
      b = dotprod(a2, ac[i]);
      c = dotprod(a3, ac[i]);
      ac[i].x = -b;
      ac[i].y = -a;
      ac[i].z = -c;
   }
   
} /* End of convert */

/**********************************************************************
**********************************************************************/
void second(Vector3 ai[], double tau, int natom, int n1, int n2, 
            int n3, int n4, Vector3 ao[])

{  int       i, j;
   double    del1, del2, d, snk, csk, dc, snw, csw;
   Vector3   bb;
   Matrix    s, t;

   tau *= DEGtoRAD;
   del1 = delta(ai[n1], ai[n2], ai[n4]);
   del2 = delta(ai[n4], ai[n3], ai[n1]);
   d = dist(ai[n1], ai[n4]);
   del1 *= DEGtoRAD;
   del2 *= DEGtoRAD;
   tau = M_PI + del2 - (tau + del1);
   snk = sin(tau);
   csk = cos(tau);
   dc = dist(ai[n3], ai[n4]);
   dc = 1.0 / dc;
   snw = (ai[n3].z * dc) / sin(del2);
   csw = sqrt(1.0 - snw * snw);
   if (ai[n3].x < 0.0)
     csw = -csw;
   initmat(s, 0.0);
   s[0][0] = csw * csk;
   s[0][1] = csw * snk;
   s[0][2] = -snw;
   s[1][0] = -snk;
   s[1][1] = csk;
   s[1][2] = 0.0;
   s[2][0] = snw * csk;
   s[2][1] = snw * snk;
   s[2][2] = csw;
   for (i=0; i<4; i++) 
     for (j=0; j<4; j++) 
       t[i][j] = s[j][i];
   for (i=0; i<natom; i++) 
   {
     bb = ai[i];
     bb.y += d;
     TransfV(bb, t, &ao[i]);
   }

} /* End of second */

/**********************************************************************
   Ex.  1->4 linkage:

   0  ->   O(1)  ->  norg
   1  ->   C(1)  ->  nxy
   2  ->   C(4)  ->  nxyz
   3  ->   O(4)  ->  ny
   4  ->   O(5)  ->  nh
   5  ->   C(5)2 ->  nhpr
**********************************************************************/
void eneh(int norg, int nxy, int nxyz, int ny, int nh, int nhpr,
          double ang1, double ang2, double tau, 
          double *na, double *n, double *h, double *doo)

{  int     i, j, k;
   double  del1, del2, d, snw, csw, snk, csk;
   double  dphi, dpsi, phi, phi1, anf1, anf2, ffi, fi;
   double  csth, snth, delt;
   Vector3 bx[6], ax[6], cx[6], ex, ee, ek, eh;
   Matrix  r1, r, a, q;

   bx[0] = A[norg].pos;
   bx[1] = A[nxy].pos;
   bx[2] = A[nxyz].pos;
   bx[3] = A[ny].pos;
   bx[4] = A[nh].pos;
   bx[5] = A[nhpr].pos;

/* for (i=0; i<6; i++)
     printf(" bx[%d]: %10.5lf %10.5lf %10.5lf\n", i, bx[i].x, bx[i].y, bx[i].z);
*/
   convert(bx, 6, 0, 3, 1, ax);
/* for (i=0; i<6; i++)
     printf(" ax[%d]: %10.5lf %10.5lf %10.5lf\n", i, ax[i].x, ax[i].y, ax[i].z);
*/
//// 20030118 pih: &cx -> cx [compile is clean, no warnings]
   second(ax, tau, 6, 0, 1, 2, 3, cx);
/* for (i=0; i<6; i++)
     printf(" cx[%d]: %10.5lf %10.5lf %10.5lf\n", i, cx[i].x, cx[i].y, cx[i].z);
*/
   dphi = torval(ax[4], ax[1], ax[0], cx[2]);
   dpsi = torval(ax[1], ax[0], cx[2], cx[5]);
   dpsi += 180.0;

   anf1 = ang1 - dphi;
   anf2 = ang2 - dpsi;

   del1 = delta(ax[0], ax[1], ax[3]);
   del2 = delta(ax[3], ax[2], ax[0]);
   /* printf("DEL1: %10.2lf   DEL2: %10.2lf\n", del1, del2); */
   uvec(ax[1], ax[0], &ee);
   d = dist(ax[2], ax[3]);
   del2 *= DEGtoRAD;
   snw = ax[2].z / ( d * sin(del2));
   d = dist(ax[0], ax[3]);
   csw = sqrt(1.0 - snw*snw);
   if (ax[2].x < 0.0)
      csw = -csw;
   del1 *= DEGtoRAD;
   uvec(ax[2], ax[3], &ek);
   ek.x *= -1.0;
   ek.y *= -1.0;
   ek.z *= -1.0;
   tau *= DEGtoRAD;
   tau = M_PI + del2 - (tau + del1);
   snk = sin(tau);
   csk = cos(tau);
   initmat(r1, 0.0);
   initmat(r, 0.0);
   initmat(q, 0.0);
   initmat(a, 0.0);
   r1[0][0] = csw * csk;
   r1[0][1] = csw * snk;
   r1[0][2] = -snw;
   r1[1][0] = -snk;
   r1[1][1] = csk;
   r1[1][2] = 0.0;
   r1[2][0] = snw * csk;
   r1[2][1] = snw * snk;
   r1[2][2] = csw;

   phi = anf1 * DEGtoRAD * 0.5;
   mxrot(ee, phi, r);
   for (i=0; i<3; i++) 
     for (j=0; j<3; j++) 
     {  q[i][j] = 0.0;
        for (k=0; k<3; k++)
          q[i][j] += r1[i][k]*r[k][j];
     }
   ffi = ang2;
   fi  = anf2;
   phi1 = (fi + 180.0)*DEGtoRAD*0.5;
   mxrot(ek, phi1, r);  
   for (i=0; i<3; i++) 
     for (j=0; j<3; j++) 
     {  a[i][j] = 0.0;
        for (k=0; k<3; k++)
          a[i][j] += r[i][k]*q[k][j];
     }
   csth = 0.5 * (a[0][0] + a[1][1] + a[2][2] - 1.0);
   /* printf("csth = %.3lf\n", csth); */
   if (fabs(csth) > 1.0) 
     if (csth > 0.0)
       csth = 1.0;
     else
       csth = -1.0;
   delt = acos(csth);
   snth = -1.0 / (2.0 * sin(delt));
   *na= delt * RADtoDEG;
   *n = 360.0 / (delt * RADtoDEG);
   eh.x = ( a[2][1] - a[1][2] ) * snth;
   eh.y = ( a[0][2] - a[2][0] ) * snth;
   eh.z = ( a[1][0] - a[0][1] ) * snth;
   *h = eh.y*d;
   *doo = d;
   nfun++;
   
} /* End of eneh */

/**********************************************************************
   Some global static variables needed to make calc_dn2 a function 
   of only p[]
**********************************************************************/
#define MAXDF  12
static  int    maxrtv, RTI[MAXDF];
static  int    A1, A2, A3, A4, A5, A6;
static  double no;

/**********************************************************************
	Calculate the helical parameters as a function of RTV[]

	This function will return the value of:

	no + (rpt - no)*(rpt - no)
**********************************************************************/
double calc_dn2(double *p)

{  int     i;  
   double  tau;
   double  na, rpt, apr, doo, dn2;
   Tors    Pi;
 
   Pi = RT[RTI[0]];
   tau = angval(A[Pi.over1].pos,A[Pi.over2].pos,A[Pi.to].pos);

   for (i=0; i<maxrtv; i++)
   {  Pi = RT[RTI[i]];
      rotbond(Pi, p[i+1]);
      if (print) printf("%10.3lf", p[i+1]);
   }
   if (print) printf(" || ");

   eneh(A1,A2,A3,A4,A5,A6,p[1],p[2],tau,&na,&rpt,&apr,&doo);

   if (print) printf("%10.3lf  %10.3lf  %10.7lf\n", doo, apr, rpt); 

   dn2 = (rpt - no)*(rpt - no);
   return(no+dn2);

} /* End of calc_dn2 */

/**********************************************************************
**********************************************************************/
void optimize_n(int phino, int psino, int omeno, int fres, int lres) 

{  int       i, j;
   int       iter, res1, res2;
   double    delta = 0.10, tole = 1.0E-10;
   double    n, dn2, fret, f1, *p, *g;
   Tors      Pi;


   printf("Orthogonal Optimization of Helical Structure:\n\n");
   printf("Search for nearest iso-(n=%.5lf) contour\n",no);
   /* insert PHI, PSI and OME in list of main chain rotable bonds */
   j = 0;
   RTI[j++] = phino;
   RTI[j++] = psino;
   if (omeno !=-1)
     RTI[j++] = omeno;

   /* Create list of additional main chain rotable bonds */
   for (i=0; i<M.nrot; i++)
     if ((i != phino) AND (i != psino) AND (i != omeno))
     {  Pi = RT[i];
        res1 = A[Pi.from].res;
        res2 = A[Pi.to].res;
        if (res1 != res2)
          if ((res1 >= fres) AND (res1 <= lres))
            if ((res2 >= fres) AND (res2 <= lres))
              if ((Res[res1].status == 0) AND (Res[res2].status == 0))
                RTI[j++] = i;
   }
   maxrtv = j;

   /* create some vectors */
   p=dvector(1,maxrtv);
   g=dvector(1,maxrtv);

   printf("Relaxed main chain dihedrals:\n\n");
   for (j=0; j<maxrtv; j++)
   {
     Pi = RT[RTI[j]];
     p[j+1]=torval(A[Pi.from].pos,A[Pi.over1].pos,A[Pi.over2].pos,A[Pi.to].pos);
     printf("\tROT[%2d]: %s[%2d] - %s[%2d] --> %8.3lf\n", RTI[j], 
             A[Pi.over1].lab, Pi.over1, A[Pi.over2].lab, Pi.over2, p[j+1]);
   }
   printf("\n");

   /* calculate the energy and its derivatives */
   fret = calc_dn2(p);
   for (i=1; i<=maxrtv; i++)
   {  p[i] += delta;
      f1 = calc_dn2(p);
      p[i] -= delta;
      g[i] = -(f1 - fret)/delta;
   }

   /* perform line search to find nearest orthogonal iso-N */
   for (i=1; i<=maxrtv; i++)
     printf("----%2d----",i);
   printf(" || ");
   printf("-----doo----------H-------N-------\n"); 

   print = TRUE;
   linmin(p,g,maxrtv,EPSILON,&fret,calc_dn2);
   print = FALSE;

   printf("\nRefined main chain vector:\n\n", fret);
   for (j=0; j<maxrtv; j++)
     printf("\tROT[%2d]: %s[%2d] - %s[%2d] --> %15.7lf\n", RTI[j], 
             A[RT[RTI[j]].over1].lab, RT[RTI[j]].over1, 
             A[RT[RTI[j]].over2].lab, RT[RTI[j]].over2, p[j+1]);

   /* free heap memory for vectors */
   free_dvector(p,1,maxrtv);
   free_dvector(g,1,maxrtv);

} /* End of optimize_n */

/**********************************************************************
   This rutine will setup the calculations of helical parameters.

   Three repeat units should have been build:

   [ -A-B-C- ] 3
        |
        B'

   The polymer repeat fragment will now be treated as a virtual fragment:

        _ _ _ _ _
       /         \
 [ -O-C           C-O- ]n
       \_ _ _ _ _/

     g2nr       g2r
**********************************************************************/
void calc_helix(BOOLEAN optimize, double nfold)

{  int      i;
   int      rpr, sub_from, sub_to;
   int      g2r, g2nr, offset=0;
   int      phino=-1, psino=-1, omeno=-1;
   double   na=0.0, rpt=0.0, apr=0.0, apt=0.0, doo=0.0;
   double   phi, psi, ome, tau;
   Matrix   MM;
   Tors     Pi;

   printf("Calculation of Helical Parameters\n\n");
   printf("\tNote: 3 repeats must have been build\n");
   printf("\tNote: A table of rotatable bonds must have been created\n\n");

   /* Check if list of rotatable bonds has been created */
   if (M.nrot == 0)
   {  printf("A list of rotatable bonds have not been created\n");
      printf("No helical calculations will be done !!!\n\n\n");
      exit(0);
   }
   /* Get number of residues per repeat sequence */
   if (((M.nres - 1) % 3) != 0)
   {  printf("Number of builded repeats not equal to 3\n");
      printf("No helical calculations will be done !!!\n\n\n");
      exit(0);
   }
   rpr = (M.nres-1) / 3;

   /*   Find atomic indexes                */
   /*   PHI:           A5 - A2 - A1 -(A3)  */
   /*   PSI:          (A2)- A4 - A3 - A6   */
   /*   OME:           A4 - A3 - A6 -(A5)  */
   /*   TAU:           A2 - A1 -(A3)       */

   /* find numbers of non-reducing and reducing residues on G2 */
   /* BEWARE: sidechains can mess up this definition */
   g2nr = rpr + 1;
   g2r  = rpr + rpr;
   if (Res[g2r].status != 0)
   {  for (i=g2r; i>g2nr; i--)
        if (Res[i].status == 0)
          break;
      offset = g2r - i;
      g2r = i;
   }
   printf("Main chain repeat segment: %d --> %d  (offset = %d)\n",
             g2nr, g2r, offset);
   
   /* get linkage definition between repeats */
   sub_from = Res[g2r].from;
   sub_to   = Res[g2r].to;
   phi      = Res[g2r].phi;
   psi      = Res[g2r].psi;
   tau      = Res[g2r].theta;

   /* reducing end */
   switch(sub_from)
   {  case 1: if (Res[g2r].ringtype==5)         /* aldopyranose */
              {  A5 = AtomIndex(g2r, "O_4");
                 A2 = AtomIndex(g2r, "C_1");
                 A1 = AtomIndex(g2r, "O_1");
              }
              else
              {  A5 = AtomIndex(g2r, "O_5");   /* pyranose */
                 A2 = AtomIndex(g2r, "C_1");
                 A1 = AtomIndex(g2r, "O_1");
              }
              break;
      case 2: A5 = AtomIndex(g2r, "O_6");
              A2 = AtomIndex(g2r, "C_2");
              A1 = AtomIndex(g2r, "O_2");
              break;
     default: printf("\nHelix: REDUCING ERROR\n");
              exit(0);
   }

   /* nonreducing end */
   switch(sub_to)
   {  case 2: A4 = AtomIndex(g2nr-offset-1, "O_1");
              A3 = AtomIndex(g2nr, "C_2");
              A6 = AtomIndex(g2nr, "C_3");
              break;
      case 3: A4 = AtomIndex(g2nr-offset-1, "O_1");
              A3 = AtomIndex(g2nr, "C_3");
              A6 = AtomIndex(g2nr, "C_4");
              break;
      case 4: A4 = AtomIndex(g2nr-offset-1, "O_1");
              A3 = AtomIndex(g2nr, "C_4");
              A6 = AtomIndex(g2nr, "C_5");
              break;
      case 5: /* aldofuranose only */
              A4 = AtomIndex(g2nr-offset-1, "O_1");
              A3 = AtomIndex(g2nr, "C_5");
              A6 = AtomIndex(g2nr, "C_4");
              break;
      case 6: A4 = AtomIndex(g2nr-offset-1, "O_1");
              A3 = AtomIndex(g2nr, "C_6");
              A6 = AtomIndex(g2nr, "C_5");
              break;
      case 8: A4 = AtomIndex(g2nr-offset-1, "O_1");
              A3 = AtomIndex(g2nr, "C_8");
              A6 = AtomIndex(g2nr, "C_7");
              break;
     default: printf("\nHelix: NONREDUCING ERROR\n");
              exit(0);
   }
   if (debug)
   {  printf("PHI:  %d - %d - %d -(%d)\n", A5, A2, A1, A3);
      printf("PSI: (%d)- %d - %d - %d \n", A2, A4, A3, A6);
      if (sub_to == 6)
         printf("OME:  %d - %d - %d -(%d)\n", A4, A3, A6, A5);
      printf("TAU:  %d - %d -(%d)\n", A2, A1, A3);
   }

   for (i=0; i<M.nrot; i++)
   {  Pi = RT[i];

      /* identify PHI in RT[] */
      if ((Pi.from == A5) AND (Pi.over1 == A2) AND (Pi.over2 == A1))
         phino = i;

      /* identify PSI in RT[] */
      if ((Pi.over1 == A4) AND (Pi.over2 == A3) AND (Pi.to == A6))
         psino = i;

      /* identify OME in RT[] */
      if ((Pi.over1 == A6) AND (Pi.over2 == A3) AND (Pi.to == A4))
         omeno = i;
   }
   phi = torval(A[RT[phino].from].pos,A[RT[phino].over1].pos,
                A[RT[phino].over2].pos,A[RT[phino].to].pos);
   psi = torval(A[RT[psino].from].pos,A[RT[psino].over1].pos,
                A[RT[psino].over2].pos,A[RT[psino].to].pos);
   tau = angval(A[RT[phino].over1].pos,A[RT[phino].over2].pos,
                A[RT[phino].to].pos);
   printf("\nInput for helical calculations:\n");
   printf("\tPHI[%2d]: %s - %s - %s - %s  == %7.2lf\n", 
           phino, A[A5].lab, A[A2].lab, A[A1].lab, A[RT[phino].to].lab, phi);
   printf("\tPSI[%2d]: %s - %s - %s - %s  == %7.2lf\n", 
           psino, A[RT[psino].from].lab, A[A4].lab, A[A3].lab, A[A6].lab, psi);
   printf("\tTAU:     %s - %s - %s == %7.2lf\n", 
                  A[A2].lab, A[A1].lab, A[RT[phino].to].lab, tau); 

   eneh(A1,A2,A3,A4,A5,A6,phi,psi,tau,&na,&rpt,&apr,&doo);

   apt=rpt*apr/rpr;
   printf("\n---HELICAL-PARAMETERS--------------------------\n");
   printf("Residues per repeat                : %4d\n", rpr);
   printf("Residues per helical turn       (n): %7.2lf \n", rpt);
   printf("Helical rotation per repeat    (na): %7.2lf \n", na);
   printf("----------------------------------------------\n");
   printf("Helical repeat extension  (O- - -O): %7.2lf A\n", doo);
   printf("Helical repeat advancement      (h): %7.2lf A\n", apr);
   printf("Helical fiber repeat          (n*h): %7.2lf A\n", apt);
   printf("----------------------------------------------\n\n");

   /* now optimize the chain dihedrals to give an integer n */
   /* we will minimize the number: fabs(nx - int(no))  */
   if (optimize)
   {  if (nfold != 0.0)
        no = nfold;
      else
      { if ((ceil(rpt)-rpt) < (rpt-floor(rpt)))
          no = ceil(rpt);
        else
          no = floor(rpt);
      }
      optimize_n(phino,psino,omeno,g2nr,g2r);
   }
   tau = angval(A[RT[phino].over1].pos,A[RT[phino].over2].pos,
                A[RT[phino].to].pos);
   printf("\tTAUAV:     %s - %s - %s == %7.2lf\n", 
                  A[A2].lab, A[A1].lab, A[RT[phino].to].lab, tau); 

} /* End of calc_helix */

/* End of file */
