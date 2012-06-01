/******************************************************************POLYS
** coord.c                                                            **
**                                                                    **
** Molecular Topology                                                 **
** ------------------                                                 **
** This file contains the source code for performing coordinate       **
** manipulations in cartesian space                                   **
**                                                                    **
** Written by Soren Balling Engelsen, INRA-93/94                      **
***********************************************************************/

#include <stdio.h>
#include <math.h>
#include "polys.h"
#include "extern.h"

int
getAno(int resnum, Label lab)
{
   int i;
   
   for (i=0; i<M.nat; i++)
      if (A[i].res == resnum)
         if (strcmp(lab, A[i].lab) == 0)
            return (i);
    
   return (-1);
}

void
Coord_copy()
{
   int i;

   for (i=0; (i<M.nat) AND (i<MAXATOMR); i++)
      Ar[i]  = A[i].pos;
}

void
symneighbor(int nat, Matrix TM)
{
   int i;
   Vector3 afrac;

   printf("\nAdding symmetry neighbor:\n");
   matprint(TM);
   for (i=0; i<nat; i++)
   {  A[M.nat+i] = A[i];
      TransfV(Ar[i],TM,&afrac);
      TransfV(afrac,CTF,&A[M.nat+i].pos);
   }
   M.nat += nat;  
   printf("\nNumber of atoms in main set increased to: %d\n", M.nat);
}

/*
 * atomnumber no  to origo 
 * atomnumber nx  to the X-axis
 * atomnumber nxy to the XY-plane 
 */
void
Coord_orient(no, nx, nxy)
{
   Orient(1, M.nres, no, nx, nxy);
   printf("Coordinate set has been transformed:\n\n");
   printf("At the origin:\n");
   printf("\tAtom no. %3d: %-8s [%d]\n", no, A[no].lab, A[no].res);
   printf("On the x-axis:\n");
   printf("\tAtom no. %3d: %-8s [%d]\n", nx, A[nx].lab, A[nx].res);
   printf("In the xy-plane:\n");
   printf("\tAtom no. %3d: %-8s [%d]\n", nxy, A[nxy].lab, A[nxy].res);
}

void
Coord_rms(BOOLEAN heavy)
{
   int    i, count, bigi, smali;
   double dx, dy, dz, d, ds, dmin, dmax, dxs, dys, dzs, rms;

   dxs = 0.0;
   dys = 0.0;
   dzs = 0.0;
   ds  = 0.0;
   dmin = 1000.0;
   dmax =    0.0;
   count = 0;
   for (i=0; (i<M.nat) AND (i<MAXATOMR); i++)
   {  if ((!heavy) OR ((A[i].ztyp / 100) != 1))
      {  count++;
         dx  = Ar[i].x - A[i].pos.x;
         dxs += dx*dx;
         dy  = Ar[i].y - A[i].pos.y;
         dys += dy*dy;
         dz  = Ar[i].z - A[i].pos.z;
         dzs += dz*dz;
         d  = sqrt(dx*dx + dy*dy + dz*dz);
         ds += d;
         if (d < dmin)
         {  dmin = d;
            smali = i;
         }
         if (d > dmax)
         {  dmax = d;
            bigi = i;
         }
      }
   }
   rms = sqrt((dxs + dys + dzs)/(3*count));
   dxs /= count;
   dys /= count;
   dzs /= count;
   ds  /= count;
   dxs = sqrt(dxs);
   dys = sqrt(dys);
   dzs = sqrt(dzs);

   printf("Number of atoms compared:   %d\n\n", count);
   printf("The smallest distance found was:   %.4lf %-8s\n",dmin,A[smali].lab);
   printf("The largest distance found was:    %.4lf %-8s\n",dmax,A[bigi].lab);
   printf("The average distance was:          %.4lf\n", ds);
   printf("The rms difference in x-direction: %.4lf\n", dxs);
   printf("The rms difference in y-direction: %.4lf\n", dys);
   printf("The rms difference in z-direction: %.4lf\n", dzs);
   printf("\nThe rms difference:                %.4lf\n", rms);
}

void
jacoby(Matrix A, Matrix B)
{
   int     i, j, maxi, maxj, iter, maxiter = 12;
   double  emax, ts, del, y, r, a, b, d, f1, fcos, fsin;

   initmat(B, 0.0);
   for (i=0; i<3; i++)
      B[i][i] = 1.0;
   iter = 0;
 
   do 
   {  iter++;
      /* find largest offdiagonal element in A */
      emax = 0.0;
      for (i=1; i<3; i++)
         for (j=1; j<3; j++)
            if ((i != j) AND (fabs(A[i][j]) > emax))
            {  emax = fabs(A[i][j]);
               maxi = i;
               maxj = j;
            }

      /* calculate rotation */
      ts  = A[maxi][maxj]*A[maxi][maxj];
      del = A[maxi][maxi] - A[maxj][maxj];
      r   = sqrt(4.0*ts + del*del);
      a   = sqrt(fabs((r+del)/(2.0*r)));
      if (a < 0.707)
      {  b = -a;
         a = sqrt(1.0 - b*b);
      }
      else
         b = -sqrt(1.0 - a*a);
      if (del/A[maxi][maxj] < 0)
         b = -b;

      /* Orthogonal transformation of matrix */
      for (j=0; j<3; j++)
         if (maxj != j)
            if (maxi != j)
            {  A[maxi][j] = a*A[j][maxi] - b*A[j][maxj];
               A[maxj][j] = b*A[j][maxi] + a*A[j][maxj];
            }
      for (j=0; j<3; j++)
      {  A[j][maxi] = A[maxi][j];
         A[j][maxj] = A[maxj][j];
      }
      d = A[maxi][maxi] + A[maxj][maxj];
      A[maxi][maxi] = a*a*A[maxi][maxi] + b*b*A[maxj][maxj] 
                      - 2.0*a*b*A[maxi][maxj];
      A[maxj][maxj] = d - A[maxi][maxi];
      A[maxj][maxi] = 0.0;
      A[maxi][maxj] = 0.0;

      /* enter rotation in eigenvector matrix */
      for (j=0; j<3; j++)
      {  f1 = B[j][maxi];
      }  

   } while (iter < 12);

/*

  140 ITER=ITER+1   
  310 CONTINUE       
      D=X(K,K)+X(L,L)                 
      X(K,K)=A*A*X(K,K)+B*B*X(L,L)-2.*A*B*X(K,L)
      X(L,L)=D-X(K,K)                  
      X(L,K)=0.                         
      X(K,L)=0.                          
C --- ENTER ROTATION IN EIGENVECTOR MATRIX
  320 DO 330 J=1,N        
      F1=E(J,K)            
      E(J,K) = A*F1-B*E(J,L)
  330 E(J,L) = B*F1+A*E(J,L)
      IF(ITER-NSQP) 140,390,390
  340 BIGR=0.00004
      NMI=N-1    
      DO 380 I=1,NMI
      IPI=I+1      
      DO 380 J=IPI,N
      IF(X(I,J))350,380,350
  350 DEL=X(I,I)-X(J,J)
      IF(DEL) 360,380,360
  360 ROT=ABS(X(I,J)/DEL)
      IF(ROT-BIGR)380,380,370
  370 BIGR=ROT
      K=I
      L=J
  380 CONTINUE
      IF(BIGR-0.00005) 390,390,200
  390 RETURN
      END
*/
} /* End of jacoby */

void
Coord_stat()
{
   int    i;
   double xom, yom, zom, omw, am;
   double sxx, syy, szz, sxy, syz, szx, xi, yi, zi, xx, yy, zz, abc;
   Matrix INERTI, DIAG;

   /* first calculate center of mass */
   xom = 0.0;
   yom = 0.0;
   zom = 0.0;
   omw  = 0.0;
   for (i=0; i<M.nat; i++)
   {  am = atomass(i);
      xom += A[i].pos.x*am;
      yom += A[i].pos.y*am;
      zom += A[i].pos.z*am;
      omw += am;
   }
   xom /= omw;
   yom /= omw;
   zom /= omw;
   printf("Center of mass was:\n");
   printf("   x = %.2lf   y = %.2lf   z = %.2lf\n\n", xom, yom, zom);

   /* then move origin to the center of mass */
   for (i=0; i<M.nat; i++)
   {  A[i].pos.x -= xom;
      A[i].pos.y -= yom;
      A[i].pos.z -= zom;
   }

   /* calculate moments of inertia */
   sxx = 0.0;
   syy = 0.0;
   szz = 0.0;
   sxy = 0.0;
   szx = 0.0;
   syz = 0.0;
   xx = 0.0;
   yy = 0.0;
   zz = 0.0;
   for (i=0; i<M.nat; i++)
   {  am = atomass(i);
      xi = A[i].pos.x;
      yi = A[i].pos.y;
      zi = A[i].pos.z;

      /* summing the M*X*X elements */
      xx += am*xi*xi;
      yy += am*yi*yi;
      zz += am*zi*zi;

      /* summing the M*X*Y elements */
      sxy += am*xi*yi;
      szx += am*xi*zi;
      syz += am*yi*zi;

   }
   /* calculation of the sums of M*X*X + M*Y*Y elements */
   sxx = yy + zz;
   syy = xx + zz;
   szz = xx + yy;

   /* calculation of the three principal moments of inertia */
   abc = sxx*syy*szz-2.0*syz*szx*sxy-syy*szx*szx-sxx*syz*syz-szz*sxy*sxy;

   /* building the inertiatensor */
   initmat(INERTI, 0.0);
   INERTI[0][0] =  sxx;
   INERTI[1][1] =  syy;
   INERTI[2][2] =  szz;
   INERTI[0][1] = -sxy;
   INERTI[1][0] = -sxy;
   INERTI[0][2] = -szx;
   INERTI[2][0] = -szx;
   INERTI[1][2] = -syz;
   INERTI[2][1] = -syz;
   INERTI[0][3] = 0.0;
   INERTI[1][3] = 0.0;
   INERTI[2][3] = 0.0;
   INERTI[3][3] = 1.0;

/*   matprint(INERTI); */

   /* diagonalize matrix */
   jacoby(INERTI, DIAG); 

/*   matprint(DIAG);   */

} /* End of Coord_stat */
