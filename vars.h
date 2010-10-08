/*****************************************************************************
   "vars.h"

                GLOBAL VARIABLE DECLARATIONS 

SBE93-INRA
*****************************************************************************/
BOOLEAN  debug = FALSE;  /* Boolean control variables */
BOOLEAN  plot = FALSE;
BOOLEAN  print = FALSE;
BOOLEAN  silent = FALSE;
int      lineno = 1;
int      seed = -7777;
int      maxit  = 2000;
int      maxfun = 10000;
int      nfun = 0;
int      da = 0;               /* number of atoms to be deleted         */
int      nmat;                 /* number of symmetry operations         */
int      format = 6;           /* default xyz format: MONO              */
                               /*          AMBER = 0                    */
                               /*          CHARMM= 1                    */
                               /*          CFF   = 2                    */
                               /*          CSD   = 3                    */
                               /*          DAT   = 4                    */
                               /*          MM3   = 5                    */
                               /*          MONO  = 6                    */
                               /*          PDB   = 7                    */
                               /*          PF    = 8                    */
                               /*          PLOT  = 9                    */
                               /*          SYBYL = 10                   */
                               /*          SPF   = 11                   */
                               /*          PIM   = 12                   */
                               /*          XYZ   = 13                   */
int      genopt = 0;           /* default torsion generation: ONE       */
                               /*          ONE   = 0                    */
                               /*          ALL   = 1                    */
int      space  = 0;           /* default minimization space: CARTESIAN */
                               /*          CARTE = 0                    */
                               /*          DIHED = 1                    */
double   sc14fc, sc14fv;
double   grdtol = 0.001;
double   enetol = 0.00001;
double   ttot,tbuild,tint,tpara,tnb,tbo,tv,tt,trot,tother,tid,ttest;

char *ProgName;            /* name of the program */
char *UserName;            /* name of the user */

/* Constructed datatypes */
/* orthogonal transformations */
Crystal    CELL = {1.0, 1.0, 1.0, 90.0, 90.0, 90.0};
Matrix     CTF;
Matrix     MAT[MAXIMG];      /* Symmetry operations */
Molecule   M;                /* molecular properties */
Residue    Res[MAXRES];      /* table of the residues */
atom       A[MAXATOM];       /* table of the atoms */
Vector3    Ar[MAXATOMR];     /* table of reference atoms */
Bond       B[MAXBOND];       /* table of the covalent bonds */
Angl       T[MAXTHETA];      /* table of the valence angles */
Tors       P[MAXPHI];        /* table of the torsional angles */
Tors       RT[MAXROT];       /* table of rotatable dihedrals */
Oops       O[MAXOOP];        /* table of out-of-plane angles */
Vector3    F[MAXATOM];       /* the force vector */
int        CONN[MAXATOM];    /* connection level vector */
int        INB[MAXATOM];     /* index vector for nonbonded exclusions */
int        I14[MAXATOM];     /* index vector for 1-4 nonbonded interactions */
int        ANB[3*MAXATOM];   /* nonbonded exclusions */
int        A14[3*MAXATOM];   /* nonbonded 1-4 interactions */
int        DA[MAXATOM];      /* hydrolysis vector, see 'da' */

Nbap       NBP[MAXTYPE];

/* Help message */

static char *help_message[] = {
"where options include:",
"    -display host:dpy            X server to use",
"    -geometry geom               geometry of window to use",
"    -d                           include full debug information",
"    -t                           include full timing information",
"    -np                          exclude graphical part",
"    -ran                         randomize torsional angles in linkages",
"    -prob                        probable torsional angles in linkages",
"    -helix                       helix torsional angles in linkages",
"    -u username                  name of user",
"    -m molname                   associated name of molecule/session",
NULL};


/* End of file */
