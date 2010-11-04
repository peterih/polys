/*
 * This module contains rutines for allocating vectors and matrices.
 * ADAPTED FROM NUMERICAL RECIPES
 * Written by Soren Balling Engelsen, INRA-93.
 */
#include <stdio.h>
#include <malloc.h>

static void
alerror(char error_text[])
{
   void exit();

   fprintf(stderr,"Allocation run-time error...\n");
   fprintf(stderr,"%s\n",error_text);
   fprintf(stderr,"...now exiting to system...\n");
   exit(1);
}

/*
 * Allocate a vector of integers on the heap.
 * BEWARE: the vector will start with index 1.
 */
int
*ivector(int nl, int nh)
{
   int *v;

   v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
   if (!v) alerror("allocation failure in ivector()");
   return v-nl;
}

/*
 * Allocate a vector of float on the heap.
 * BEWARE: the vector will start with index 1.
 */
float
*vector(int nl, int nh)
{
   float *v;

   v = (float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
   if (!v) alerror("allocation failure in vector()");
   return v-nl;
}

/*
 * Allocate a vector of double on the heap.
 * BEWARE: the vector will start with index 1.
 */
double
*dvector(int nl, int nh)
{
   double *v;

   v = (double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
   if (!v) alerror("allocation failure in dvector()");
   return v-nl;
}

/*
 * Allocate a matrix of integer on the heap.
 * BEWARE: the matrix will start with index [1,1].
 */
int
**imatrix(int nrl, int nrh, int ncl, int nch)
{
	int i,**m;

	m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) alerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) alerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}

/*
 * Allocate a matrix of float on the heap.
 * BEWARE: the matrix will start with index [1,1].
 */
float
**matrix(int nrl, int nrh, int ncl, int nch)
{
   int i;
   float **m;

	m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
	if (!m) alerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
		if (!m[i]) alerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

/*
 * Allocate a matrix of double on the heap.
 * BEWARE: the matrix will start with index [1,1].
 */
double
**dmatrix(int nrl, int nrh, int ncl, int nch)
{
   int i;
   double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) alerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) alerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

/* Free heap space occupied by a previously allocated ivector of integers */
void free_ivector(int *v, int nl, int nh)
{
   free((char*) (v+nl));
}

/* Free heap space occupied by a previously allocated vector of double */
void
free_vector(float *v, int nl, int nh)
{
   free((char*) (v+nl));
}

/* Free heap space occupied by a previously allocated dvector of double */
void
free_dvector(double *v, int nl, int nh)
{
   free((char*) (v+nl));
}

/* Free heap space occupied by a previously allocated imatrix of integer */
void
free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
   int i;

   for(i=nrh;i>=nrl;i--) 
      free((char*) (m[i]+ncl));
   free((char*) (m+nrl));
}

/* Free heap space occupied by a previously allocated matrix of float */
void
free_matrix(float **m, int nrl, int nrh, int ncl, int nch)
{
   int i;

   for(i=nrh;i>=nrl;i--) 
      free((char*) (m[i]+ncl));
   free((char*) (m+nrl));
}

/* Free heap space occupied by a previously allocated matrix of double */
void
free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
   int i;

   for(i=nrh;i>=nrl;i--) 
      free((char*) (m[i]+ncl));
   free((char*) (m+nrl));
}
