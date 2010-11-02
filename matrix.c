/*
 * This module contains functions for matrix manipulations for three 
 * dimensional modelling.
 * Written by Soren Balling Engelsen, INRA-93.
 */
#include <stdio.h>
#include <math.h>
#include "polys.h"

/* The length of a 3D vector */
double
length(Vector3 a)
{
   return (sqrt(a.x*a.x + a.y*a.y + a.z*a.z)); 
}

/* The dot product of two 3D vectors */
double
dotprod(Vector3 a, Vector3 b)
{
   return (a.x*b.x + a.y*b.y + a.z*b.z); 
}

/* The cross product of two 3D vectors */
void crossprod(Vector3 a, Vector3 b, Vector3 *c)
{
   c->x = a.y*b.z - b.y*a.z;
   c->y = a.z*b.x - b.z*a.x;
   c->z = a.x*b.y - b.x*a.y; 
}

/* Normalize a 3D vector */
double
normalize(Vector3 *a)
{
   double l;

   l = sqrt(a->x*a->x + a->y*a->y + a->z*a->z); 
   a->x /= l;
   a->y /= l;
   a->z /= l;  
   return (l);
}

/* Initialize a four by four matrix */
void
initmat(Matrix A, double value)
{
   double *ptr = &A[0][0];
   double *endptr = &A[3][3];

   for (; ptr <= endptr; *ptr++ = value)
   ;
}

/*
 * Calculate the translation matrix 'A'. Origin will be translated by
 * the vector '(tx,ty,tz)'.
 * Input:
 *   tx, ty, tz		translation in axis directions
 * Output:
 *   *A			the translation matrix
 */
void
tran3(double tx, double ty, double tz, Matrix A)
{
   register int i;
   register int j;

   for ( i=0; i<4; i++) 
   {  for ( j=0; j<4; j++)
         A[i][j]=0.0;
      A[i][i]=1.0;
   }
   A[0][3]=-tx;
   A[1][3]=-ty;
   A[2][3]=-tz;
}

/*
 * Calculate the scaling matrix 'A' given the scaling vector 
 * '(sx,sy,sz)'. One unit on the x-axis becomes 'sx' units etc.
 * Input:
 *   sx, sy, sz		scaling factors
 * Output:
 *   *A			the scaling matrix
 */
void
scale3(double sx, double sy, double sz, Matrix A)
{
   register int i;
   register int j;

   for ( i=0; i<4; i++)
      for ( j=0; j<4; j++)
         A[i][j]=0.0;
   A[0][0]=sx;
   A[1][1]=sy;
   A[2][2]=sz;
   A[3][3]=1.0;
}

/*
 * Calculate the rotation matrix 'A'. The axis are rotated anti-
 * clockwise through an angle 'theta' radians about an axis specified by 
 * 'axis'.
 * Input:
 *   axis=0			x-axis
 *   axis=1			y-axis
 *   axis=2			z-axis
 *   theta			angle to rotate
 * Output
 *   *A			rotation matrix
 */
void
rot3(int axis, double theta, Matrix A)
{
   register int i;
   register int j;
   int a1, a2;
   double c, s;

   for ( i=0; i<4; i++)
      for ( j=0; j<4; j++)
         A[i][j] = 0.0;
   A[axis][axis] = 1.0;
   A[3][3] = 1.0;
   a1 = ((axis+1) % 3);
   a2 = ((a1+1) % 3);
   c = cos(theta);
   s = sin(theta);
   A[a1][a1] = A[a2][a2] = c;
   A[a1][a2] = s;
   A[a2][a1] = -s;
}

/*
 * Calculate the matrix product 'C' of two matrices 'A' and 'B'.
 * Input:
 *   A			multiplicand (left operand)
 *   B			multiplier (right operand)
 * Output:
 *   *C			result matrix
 */
void
mult3(Matrix A, Matrix B, Matrix C)
{
   register int i;
   register int j;
   register int k;
   double   ab;

   for ( i=0; i<4; i++)
      for ( j=0; j<4; j++)
      {  ab = 0.0;
         for ( k=0; k<4; k++)
            ab = ab + A[i][k]*B[k][j];
         C[i][j] = ab;
      }
}

/* Print out the manipulation matrix */
void
matprint(Matrix A)
{
   int i, j;

   printf("\n");
   for ( i=0; i<4; i++)
   {  for ( j=0; j<4; j++)
         printf("%10.4lf", A[i][j]);
      printf("\n");
   }
   printf("\n");
}

/* Returns the 'angle' whose tangent is 'y/x' */
double
angle(double x, double y)
{
   if (fabs(x) < EPSILON)
      if (fabs(y) < EPSILON)
         return (0.0);
      else if (y > 0.0)
              return (M_PI*0.5);
           else return (M_PI*1.5);
   else if (x < 0.0)
           return (atan(y/x)+M_PI);
        else return (atan(y/x));
}
            
/*
 * Calculates the matrix 'A' representing the rotation of axes
 * through an angle 'phi' about a general line with base 'b' and 
 * direction 'd'.
 */
void
genrot(double phi, Vector3 b, Vector3 d, Matrix A)
{
   Matrix   F, G, H, W, FI, GI, HI, S, T;
   double  beta, theta, v;

   tran3(b.x, b.y, b.z, F);
   tran3(-b.x, -b.y, -b.z, FI);
   theta = angle(d.x, d.y);
   rot3(2, theta, G);
   rot3(2, -theta, GI);
   v = sqrt(d.x*d.x + d.y*d.y);
   beta = angle(d.z, v);
   rot3(1, beta, H);
   rot3(1, -beta, HI);
   rot3(2, phi, W);
   mult3(G, F, S);
   mult3(H, S, T);
   mult3(W, T, S);
   mult3(HI, S, T);
   mult3(GI, T, S);
   mult3(FI, S, A);
}

/* */
void
TransfV(Vector3 v, Matrix TM, Vector3 *w)
{
   w->x = TM[0][0]*v.x + TM[0][1]*v.y + TM[0][2]*v.z + TM[0][3];
   w->y = TM[1][0]*v.x + TM[1][1]*v.y + TM[1][2]*v.z + TM[1][3];
   w->z = TM[2][0]*v.x + TM[2][1]*v.y + TM[2][2]*v.z + TM[2][3];
}

/* */
void
TransfAV(atom A[], int ai, int ao, Matrix TM)
{
   int      i;
   Vector3 av;

   for (i=ai; i<ao; i++)
   {   av = A[i].pos;
       TransfV(av, TM, &A[i].pos);
   };
}

/* */
void InvMat3x3(Matrix A)
{
   int      i, i1, i2, j, j1, j2;
   double   det, adj;
   Matrix   AINV;
 
   /* calculate determinant */
   det =  A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
        + A[0][1]*(A[1][2]*A[2][0]-A[1][0]*A[2][2])
        + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);

   if (fabs(det) < EPSILON)
   {  printf("Attempted to invert singular Matrix\n"); 
      return; 
   }

   for (i=1; i<4; i++)
   {  i1=(i % 3)+1;
      i2=(i1 % 3)+1;
      for (j=1; j<4;j++)
      {  j1=(j % 3)+1;
         j2=(j1 % 3)+1;
         adj = (A[i1-1][j1-1]*A[i2-1][j2-1]-A[i1-1][j2-1]*A[i2-1][j1-1]);
         AINV[j-1][i-1] = adj/det;
      }
   }

   for (i=0; i<3; i++)
      for (j=0; j<3;j++)
         A[i][j] = AINV[i][j];
}

/* */
void
InvMat(Matrix M)
{
   int      i, j, *indx;
   double   det, *col, **Mi, **Mo;

   indx = ivector(1,3);
   col = dvector(1,3);
   Mi = dmatrix(1,3,1,3);
   for ( i=1; i<=3; i++) 
      for ( j=1; j<=3; j++)
         Mi[i][j]=M[i-1][j-1];

   ludcmp(Mi, 3, indx, &det);

   for (i=1; i<=3; i++)
   {  for (j=1; j<=3; j++)
         col[j] = 0.0;
      col[i] = 1.0;
      lubksb(Mi, 3, indx, col);
      for (j=1; j<=3; j++)
         M[i-1][j-1] = col[j];
   }

   free_dmatrix(Mi,1,3,1,3);
   free_ivector(indx,1,3);
   free_dvector(col,1,3);
}
