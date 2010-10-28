/**********************************************************************
POLYS FILE "polys.h"

   Header file for all public POLYS routines.

   Written by Soren Balling Engelsen, INRA-1993,94.
**********************************************************************/
#include <stddef.h>
#include <stdio.h>
#include "defs.h"
#include "typs.h"

/*******************/
/*    action.c     */
/*******************/
int actiontype(atom a, atom b, double *ab);

/*******************/
/*    bondtab.c    */
/*******************/
int Valency(atom A);
double dist(Vector3 a, Vector3 b);
void gen_bondtab();

/*******************/
/*    brent.c      */
/*******************/
double brent(double ax, double bx, double cx, double (*f)(double),
             double tol, double *xmin);

/*******************/
/*    calculat.c   */
/*******************/
void calc_dist(int i, int j);
void calc_angle(int i, int j, int k);
void calc_tors(int i, int j, int k, int l);
void dist_tabel();
double atomass(int i);
void calc_mass();
void calc_rgyr();
void calc_vol();

/*******************/
/*    coord.c      */
/*******************/
void symneighbor(int nat, Matrix TM);
void Coord_orient(int no, int nx, int nxy);
void Coord_rms(BOOLEAN heavy);
void Coord_stat();

/*******************/
/*    driver.c     */
/*******************/
int FindConn(int na, int nb, int ai, int ae, int *rotlevel);
void rotbond(Tors Pi, double torang);
void phipsidriv(int phinr, double sphi, double ephi, double incphi,
                int psinr, double spsi, double epsi, double incpsi,
                BOOLEAN silent, int option, char *fname);

/*******************/
/*    expand.c     */
/*******************/
void addmacro(int polyres, int con_b, int con_f, int macrores, char *mf,
              double phi, double psi, double omega);
void tertiary(int resmax);
void modify_val(int res, int back, int forward, double theta, 
                double phi, double psi, double omega);
void expand(int newres, char *fname, int back, int forward, 
            double theta, double phi, double psi, double omega, 
            int conres, int status);
void modify_con(int res, int conres, int status);
void copy_res(int source, int destin, int offset);
void res_list(int maxres);

/*******************/
/*    f1dim.c      */
/*******************/
double f1dim(double x);

/*******************/
/*    fileoper.c   */
/*******************/
extern FILE *fileopen(char *name, char *mode);
extern void fileclose(char *name, FILE *fp);

/*******************/
/*    frac2carte.c */
/*******************/
void frac2carte(Crystal UC, BOOLEAN ref);

/*******************/
/*    generate.c   */
/*******************/
void Generate(int option);

/*******************/
/*    initiate.c   */
/*******************/
void Initialize(void);

/*******************/
/*    intcoord.c   */
/*******************/
double angval(Vector3 a, Vector3 b, Vector3 c);
void gen_angletab();
double torval(Vector3 a, Vector3 b, Vector3 c, Vector3 d);
void gen_tortab();
void gen_nonbond();
void list_intcoor(int option);

/*******************/
/*    invmat.c     */
/*******************/
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);

/*******************/
/*    konjug.c     */
/*******************/
void cartemini(BOOLEAN silent, int itmax, double ftol);

/*******************/
/*    glyclink.c   */
/*******************/
int AtomIndex(int resnr, char *lab);
void Hydrolysis(int nwa);
void GlycLink(int red1, int redres, int red2,
              int nred1, int nredres, int nred2,
              int redmax, int redcon, int nredcon, 
              double phi, double psi, double omega);
void RepeatStruc(int rep, int bres, int eres);
void BuildInfo();

/*******************/
/*    linmin.c     */
/*******************/
void linmin(double p[], double xi[], int n, double tole, 
              double *fret, double (*func)());

/*******************/
/*    linsc.c      */
/*******************/
void linsch(int n, int *flag, double *eps, double *la, double *y,
             double fm[], double h[], double x[], double (*func)());

/*******************/
/*    matrix.c     */
/*******************/
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
double angle(double x, double y);
void genrot(double phi, Vector3 b, Vector3 d, Matrix A);
void TransfV(Vector3 v, Matrix TM, Vector3 *w);
void TransfAV(atom A[], int ai, int ao, Matrix TM);
void InvMat3x3(Matrix A);
void InvMat(Matrix M);

/*******************/
/*    minimize.c   */
/*******************/
void Minimize(int space, int maxfun, int maxit, double enetol,
              double grdtol, BOOLEAN silent);

/*******************/
/*    mnbrak.c     */
/*******************/
void mnbrak(double *ax, double *bx, double *cx,
            double *fa, double *fb, double *fc, double (*func)(double));

/*******************/
/*    orient.c     */
/*******************/
void Orient(int ia, int oa, int no, int nx, int nxy);

/*******************/
/*    palloc.c     */
/*******************/
int *ivector(int nl, int nh);
float *vector(int nl, int nh);
double *dvector(int nl, int nh);
int **imatrix(int nrl, int nrh, int ncl, int nch);
float **matrix(int nrl, int nrh, int ncl, int nch);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
void free_ivector(int *v, int nl, int nh);
void free_vector(float *v, int nl, int nh);
void free_dvector(double *v, int nl, int nh);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);

/*******************/
/*    phimin.c     */
/*******************/
void gen_rotvec();
void phimin(BOOLEAN silent,int maxfun,int maxit,double enetol,double grdtol);

/*******************/
/*    polys.c      */
/*******************/

/*******************/
/*    potener.c    */
/*******************/
double potener(int option);
double RMSE(int no, Vector3 V[]);
void printener(BOOLEAN deriv);
void tester();

/*******************/
/*    qsort.c      */
/*******************/
void Qsort(int low, int high, int ia[]);

/*******************/
/*    r_charmm.c   */
/*******************/
void read_CHARMM(char *fname);

/*******************/
/*    r_dat.c      */
/*******************/
void read_DAT(char *fname);

/*******************/
/*    r_mm3.c      */
/*******************/
void read_MM3(char *fname);

/*******************/
/*    r_mono.c     */
/*******************/
int read_MB(char *fname, BOOLEAN monob);
void read_NMB(char *fname);

/*******************/
/*    r_polys.c    */
/*******************/
void read_PF(char *fname, BOOLEAN ref);

/*******************/
/*    r_sybyl.c    */
/*******************/
void read_SYBYL(char *fname);

/*******************/
/*    random.c     */
/*******************/
double get_random(int *idum);
float ran3(int *idum);

/*******************/
/*    readcoor.c   */
/*******************/
void ReadCoord(int format, char *fname, BOOLEAN reference);

/*******************/
/*    readpara.c   */
/*******************/
void ReadParam(void);

/*******************/
/*    sinput.c     */
/*******************/
/* enum itype {INTEGER, REAL, STRING};*/
enum error {NOERROR, NOINPUT, NOSPACE, BADNUMBER, BADTYPE};
union value       /* can hold a integer, a double or a pointer to a char */
{  int    i;
   double d;
   char   *s;
};

int enter_name(char *prompt, char *name);
int get_line(FILE *fp, char *startptr, int nmax);
void putline(FILE *fp, char *startptr);
BOOLEAN yesorno(char *prompt);
int getinteger(char *prompt);
double getreal(char *prompt);
char *getstring(char *prompt);
void wline(char ctyp, unsigned short width);

/*******************/
/*    timer.c      */
/*******************/
double cput(void);
void dtime(int *year, int *month, int *day, int *hour, int *min, int *sec);
void dclock(double *t);
double etime(double tid);
void ptime(double t);
void pdate(double t);

/*******************/
/*    vbond.c      */
/*******************/
double bondcode(Bond Bi);
double paircode(int it, int jt);
void read_bp(int nbond);
double Vbond(int nbond);

/*******************/
/*    vcon.c       */
/*******************/
void read_cp(int nphi, int *ncon);
double Vcon(int nphi, int ncon);

/*******************/
/*    vnonbond.c   */
/*******************/
void neutralize(int nat);
void read_nbp(int nat, int nbond);
double Vnonbond(int derivtyp, int nat, double *Vvdw, double *Vc);

/*******************/
/*    vphi.c       */
/*******************/
void read_tp(int nphi);
double Vphi(int nphi);

/*******************/
/*    vtheta.c     */
/*******************/
void read_ap(int ntheta);
double Vtheta(int ntheta);

/*******************/
/*    writcoor.c   */
/*******************/
void WritCoord(int format, char *fname, int descarte);
void write_MB(char *fname);
void write_NMB(char *fname);
void write_PF(BOOLEAN descarte, char *fname);
void write_MM3(char *fname);
void write_DAT(BOOLEAN descarte, char *fname);
void write_SYBYL(char *fname);
void write_PIM(char *fname);
void write_SPF(char *fname);
void write_XYZ(char *fname);

/* End of file */
