/*****************************************************************************
"typs.h"

                GLOBAL DATATYPE DECLARATIONS

SBE93-INRA
*****************************************************************************/
#include <stdio.h>

typedef double Matrix[4][4];

typedef char Label[8];

/* enum ftype {AMBER,CHARMM,CFF,CSD,DAT,MM3,MONO,PDB,PF,PLOT,SYBYL,SPF,PIM,XYZ};
enum gtype {ONE,ALL};
enum stype {CARTE,DIHED}; */

typedef struct
{  int   x, y;
} Pixvec;

typedef struct
{  double   x, y;
} Vector2;

typedef struct
{  double   x, y, z;
} Vector3;

typedef struct
{  double   a, b, c, alpha, beta, gamma;
} Crystal;

typedef struct
{  char     resid[20];
   int      ringtype;
   int      status;
   int      from, to, conres;
   double   theta,phi,psi,omega;
   FILE     *fp;
} Residue;

typedef struct
{  Vector3  pos;
   double   chg;
   int      ztyp, res;
   Label    lab;
} atom;

typedef struct
{  int      from, to;
   double   Kb, alpha, l0;
} Bond;

typedef struct
{  int      from, over, to;
   double   Kt, t0;
} Angl;

typedef struct
{  int      from, over1, over2, to;
   double   K1, p1, K2, p2, K3, p3;
   int      fold;
} Tors;

typedef struct
{  int      p1, p2, p3, oop;
} Oops;

typedef struct
{  char     *id;                      /* polysaccharide id */
   char     *resid;                   /* id of last residue */
   char     *titel;                   /* name, references */
   int      nat;                      /* number of atoms */
   int      nar;                      /* number of atoms in last residue */
   int      nbond;                    /* number of bonds */
   int      ntheta;                   /* number of valence angles */
   int      nphi;                     /* number of torsional angles */
   int      ncon;                     /* number of torsional constraints */
   int      nrot;                     /* number of torsional constraints */
   int      noop;                     /* number of out-of-plane angles */
   int      nres;                     /* number of residues */
/* atom     A[ATOMMAX]; */
/* Bond     B[BONDMAX]; */
/* Angl     T[THETAMAX]; */
/* Tors     P[PHIMAX]; */
/* Oops     O[OOPMAX]; */
   double   Etot, Ekin, Epot, Eb, Et, Econ, Ep, Eo, Enbs, Evdw, Ec;
   double   temp, dipm, mw;
} Molecule;

typedef struct
{  double   PA, PB, PC;
}  Nbap;

/* End of file */
