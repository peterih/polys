/*****************************************************************************
   "extern.h"

                GLOBAL VARIABLE DECLARATIONS 

SBE93-INRA
*****************************************************************************/
extern BOOLEAN debug;  /* Boolean control variables */
extern BOOLEAN plot;
extern BOOLEAN print;
extern BOOLEAN silent;
extern int lineno;
extern int seed;
extern int maxit;
extern int maxfun;
extern int nfun;
extern int nmat;
extern double grdtol;
extern double enetol;
extern double ttot,tbuild,tint,tpara,tnb,tbo,tv,tt,trot,tother,tid,ttest;
extern double sc14fv;
extern double sc14fc;
extern int da;
extern int format;
extern int genopt;
extern int space;

extern char *UserName;            /* name of the user */
extern char *ProgName;            /* name of the program */

/* Constructed datatypes */
extern Crystal  CELL;             /* Unit cell dimensions */
extern Matrix   CTF;              /* Unit cell dimensions */
extern Matrix   MAT[MAXIMG];      /* Symmetry operations */
extern Molecule M;                /* molecular properties */
extern atom     A[MAXATOM];       /* table of the atoms */
extern Residue  Res[MAXRES];      /* table of the residues */
extern Vector3  Ar[MAXATOMR];     /* table of reference atoms */
extern Bond     B[MAXBOND];       /* table of the covalent bonds */
extern Angl     T[MAXTHETA];      /* table of the valence angles */
extern Tors     P[MAXPHI];        /* table of the torsional angles */
extern Tors     RT[MAXROT];       /* table of rotatable dihedrals */
extern Oops     O[MAXOOP];        /* table of out-of-plane angles */
extern Vector3  F[MAXATOM];       /* the force vector */
extern int      CONN[MAXATOM];    /* connection level vector */
extern int      INB[MAXATOM];     /* index vector for nonbonded exclusions */
extern int      I14[MAXATOM];     /* index vector for 1-4 interactions */
extern int      ANB[3*MAXATOM];   /* nonbonded exclusions */
extern int      A14[3*MAXATOM];   /* nonbonded 1-4 interactions */
extern Nbap     NBP[MAXTYPE];     /* nonbonded parameters */
extern int      DA[MAXATOM];      /* index vector for hydrolysis */

/* End of file */
