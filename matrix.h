/*
 * functions for matrix manipulations for three dimensional modelling.
 * SBE93-INRA
 */

double length(Vector3 a);
double dotprod(Vector3 a, Vector3 b);
void crossprod(Vector3 a, Vector3 b, Vector3 *c);
double normalize(Vector3 *a);
void initmat(Matrix A, double value);
void tran3(double tx, double ty, double tz, Matrix A);
void scale3(double sx, double sy, double sz, Matrix A);
void rot3(int axis, double theta, Matrix A);
void mult3(Matrix A, Matrix B, Matrix C);
void matprint(Matrix A);
void genrot(double phi, Vector3 b, Vector3 d, Matrix A);
void TransfV(Vector3 v, Matrix TM, Vector3 *w);
void TransfAV(atom A[], int ai, int ao, Matrix TM);
