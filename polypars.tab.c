/* A Bison parser, made from polypars
   by GNU bison 1.35.  */

#define YYBISON 1  /* Identify Bison output.  */

# define	REAL	257
# define	INTEGER	258
# define	QSTRING	259
# define	NAME	260
# define	NL	261
# define	ATTACH	262
# define	BUILD	263
# define	CALC	264
# define	COORD	265
# define	CREATE	266
# define	DATE	267
# define	DOUBLEHELIX	268
# define	END	269
# define	ENER	270
# define	GENERATE	271
# define	INIT	272
# define	HELIX	273
# define	HELP	274
# define	LABEL	275
# define	MINI	276
# define	PRIMARY	277
# define	READ	278
# define	RECORD	279
# define	MODIFY	280
# define	SET	281
# define	TIME	282
# define	WRITE	283
# define	STOP	284
# define	NSETS	285
# define	SAVE	286
# define	DIST	287
# define	DTAB	288
# define	MASS	289
# define	VOLUME	290
# define	RGYR	291
# define	COPY	292
# define	ORIENT	293
# define	RMS	294
# define	STAT	295
# define	DERIV	296
# define	NFOLD	297
# define	OPTIMIZE	298
# define	NREPS	299
# define	NSTR	300
# define	SHIFT	301
# define	ONE	302
# define	ALL	303
# define	BOND	304
# define	NBOND	305
# define	ANGLE	306
# define	TORSION	307
# define	ROTBOND	308
# define	CARTE	309
# define	DIHED	310
# define	FRAC	311
# define	GRDTOL	312
# define	ENETOL	313
# define	MAXIT	314
# define	MAXFUN	315
# define	SILENT	316
# define	NOISE	317
# define	SYMMETRY	318
# define	INTERNAL	319
# define	PARAM	320
# define	CHARGES	321
# define	REF	322
# define	PHI	323
# define	PHIB	324
# define	PHIE	325
# define	PHIINC	326
# define	PSI	327
# define	PSIB	328
# define	PSIE	329
# define	PSIINC	330
# define	RECOPT	331
# define	DEBUG	332
# define	ID	333
# define	TITLE	334
# define	RANSEED	335
# define	AMBER	336
# define	CHARMM	337
# define	CFF	338
# define	CSD	339
# define	DAT	340
# define	MM3	341
# define	MONO	342
# define	PDB	343
# define	PF	344
# define	PLOT	345
# define	SYBYL	346
# define	SPF	347
# define	PIM	348
# define	XYZ	349

#line 10 "polypars"

#include <stdio.h>
#include <string.h>
#include "defs.h"
#include "polys.h"

/*****************************************************************
   parser variables needed for the subcommands
*****************************************************************/
Label    lab;
BOOLEAN  descarte = TRUE;       /* write as cartesian coordinates */
BOOLEAN  reference = FALSE;     /* read into the reference coordinate set */
BOOLEAN  optimize = FALSE;      /* helical optimization */
int      parspace = 0;             /* default space: cartesian */
int      parformat = 6;            /* default parformat: PF */
int      i, j, no;
int      inx=0, Aijkl[4];
int      pargenopt = 1;
int      intopt = 0;
int      recopt = 0;
int      deriv = 0;
int      phino = 0;             /* no. of phi i rotb tabel */
int      psino = 0;             /* no. of psi i rotb tabel */
int      rnat = 0;              /* no. of atoms in reference set */
int      symno = 0;             /* no. of symmetry direction */
int      resno = 0;
int      nsets = 1;
int      nreps = 0;
int      nstr = 1;
double   shift = 0.0;
double   nfold = 0.0;           /* helical n parameter */
double   xmult, xtransl;
double   psib = 0.0;
double   psie = 350.0;
double   psiinc = 10.0;
double   phib = 0.0;
double   phie = 350.0;
double   phiinc = 10.0;
double   thev = 999.9;
double   phiv = 999.9;
double   psiv = 999.9;
double   omev = 999.9;
Matrix   TM;

/*****************************************************************
   parser variables needed for the BUILD command
*****************************************************************/
char     *fname;               /* the name of the MONOBANK member */
char     *mapfile;             /* the name of the DIBANK member */
char     *macrof;              /* the name of the macro file */
int      con_b, con_f;         /* glycosidic connection */
int      polyres = 0;          /* no of residue on segment to connect 
                                  sidechain */
int      repno=0;              /* current repeat sequence */
int      offset=0;             /* repeat offset */
int      status=0, oldstatus;  /* chain status */
struct  
{ int    rep;                  /* repeating factor */
  int    first;                /* first residue in repeating segment */
  int    fconnect;             /* connecting residue in repeating segment */
  int    last;                 /* last residue in repeating segment */
} reptab[10];  /* repeat table */

/*****************************************************************
   External variabels  needed for the parser
*****************************************************************/
extern   BOOLEAN debug;
extern   int lineno;
extern   int seed;
extern   int da;
extern   int silent;
extern   int maxit;
extern   int maxfun;
extern   double enetol;
extern   double grdtol;
extern   double tid;
extern   Molecule M;

#line 88 "polypars"
#ifndef YYSTYPE
typedef union { double   real;
         int      integer;
         char     *string;
         int      cmd;
       } yystype;
# define YYSTYPE yystype
# define YYSTYPE_IS_TRIVIAL 1
#endif
#ifndef YYDEBUG
# define YYDEBUG 0
#endif



#define	YYFINAL		281
#define	YYFLAG		-32768
#define	YYNTBASE	112

/* YYTRANSLATE(YYLEX) -- Bison token number corresponding to YYLEX. */
#define YYTRANSLATE(x) ((unsigned)(x) <= 349 ? yytranslate[x] : 160)

/* YYTRANSLATE[YYLEX] -- Bison token number corresponding to YYLEX. */
static const char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,    99,     2,     2,     2,     2,
     102,   104,     2,   109,   111,   110,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   103,   105,
      97,    96,    98,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   100,     2,   101,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     106,   107,   108,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     3,     4,     5,
       6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    50,    51,    52,    53,    54,    55,
      56,    57,    58,    59,    60,    61,    62,    63,    64,    65,
      66,    67,    68,    69,    70,    71,    72,    73,    74,    75,
      76,    77,    78,    79,    80,    81,    82,    83,    84,    85,
      86,    87,    88,    89,    90,    91,    92,    93,    94,    95
};

#if YYDEBUG
static const short yyprhs[] =
{
       0,     0,     1,     4,     6,     9,    15,    18,    22,    26,
      30,    33,    36,    39,    42,    45,    49,    51,    54,    56,
      59,    63,    68,    71,    73,    76,    79,    82,    89,    93,
      97,   101,   105,   108,   111,   114,   118,   122,   126,   128,
     130,   134,   138,   141,   143,   144,   147,   149,   154,   158,
     160,   166,   168,   172,   174,   176,   178,   180,   182,   184,
     190,   198,   202,   208,   210,   212,   214,   218,   220,   222,
     223,   224,   227,   230,   232,   234,   235,   238,   240,   243,
     245,   248,   250,   251,   254,   256,   257,   259,   261,   263,
     264,   267,   271,   275,   279,   280,   283,   287,   291,   293,
     294,   297,   305,   307,   309,   310,   313,   316,   319,   322,
     324,   325,   327,   330,   332,   334,   336,   339,   341,   342,
     344,   346,   348,   350,   352,   354,   356,   358,   360,   362,
     364,   366,   368,   370,   372,   374,   376,   378,   380,   382,
     384,   385,   388,   390,   394,   398,   402,   406,   408,   410,
     412,   413,   416,   418,   420,   422,   424,   426,   428,   430,
     432,   434,   436,   438,   440,   442,   444,   446,   448,   450,
     452,   454,   455,   458,   462,   466,   470,   474,   478,   482,
     486,   490,   494,   496,   498,   500,   502,   503,   506,   508,
     510,   512,   514,   516,   518
};
static const short yyrhs[] =
{
      -1,   112,   113,     0,   159,     0,   114,   159,     0,     8,
     118,   123,   131,     5,     0,     9,   127,     0,    10,    33,
     128,     0,    10,    52,   128,     0,    10,    53,   128,     0,
      10,    34,     0,    10,    35,     0,    10,    37,     0,    10,
      36,     0,    11,   132,     0,    12,    50,   128,     0,    13,
       0,    14,   137,     0,    15,     0,    16,   134,     0,    17,
      65,   136,     0,    17,    64,   141,    30,     0,    19,   139,
       0,    18,     0,    20,   149,     0,    21,     6,     0,    22,
     150,     0,    26,    53,   128,    69,    96,     3,     0,    23,
     115,    30,     0,    24,    67,     5,     0,    24,    11,   152,
       0,    24,    65,   157,     0,    24,    66,     0,    25,   155,
       0,    27,    78,     0,    27,    79,     6,     0,    27,    81,
       4,     0,    27,    80,     5,     0,    30,     0,    28,     0,
      29,    11,   152,     0,    29,    65,   157,     0,    29,    66,
       0,     6,     0,     0,   115,   116,     0,   159,     0,    97,
     117,    98,   123,     0,    97,   117,    98,     0,    97,     0,
      97,    99,   118,    98,    98,     0,   100,     0,   101,   119,
     123,     0,     6,     0,     4,     0,     4,     0,     3,     0,
       3,     0,     3,     0,   102,   125,   103,   126,   104,     0,
     102,   125,   103,   126,   105,   124,   104,     0,   120,   105,
     121,     0,   120,   105,   121,   105,   122,     0,     5,     0,
       4,     0,     4,     0,    31,    96,     4,     0,    37,     0,
      32,     0,     0,     0,   128,   129,     0,   130,   131,     0,
       6,     0,     4,     0,     0,   132,   133,     0,    38,     0,
      39,   128,     0,    40,     0,    40,    35,     0,    41,     0,
       0,   134,   135,     0,    42,     0,     0,    48,     0,    49,
       0,    54,     0,     0,   137,   138,     0,    46,    96,     4,
       0,    45,    96,     4,     0,    47,    96,     3,     0,     0,
     139,   140,     0,    43,    96,     3,     0,    43,    96,     4,
       0,    44,     0,     0,   141,   142,     0,   143,   144,   148,
     144,   148,   144,   104,     0,   159,     0,   102,     0,     0,
     145,   144,     0,   146,   106,     0,   146,   107,     0,   146,
     108,     0,   147,     0,     0,     3,     0,   109,     3,     0,
     109,     0,   110,     0,     3,     0,   109,     3,     0,   111,
       0,     0,     9,     0,    10,     0,    11,     0,    12,     0,
      13,     0,    14,     0,    15,     0,    16,     0,    17,     0,
      18,     0,    19,     0,    20,     0,    21,     0,    22,     0,
      23,     0,    24,     0,    25,     0,    26,     0,    27,     0,
      28,     0,    29,     0,     0,   150,   151,     0,    56,     0,
      58,    96,     3,     0,    59,    96,     3,     0,    60,    96,
       4,     0,    61,    96,     4,     0,    63,     0,    62,     0,
       6,     0,     0,   152,   153,     0,   154,     0,    57,     0,
      68,     0,     5,     0,     6,     0,    82,     0,    83,     0,
      84,     0,    85,     0,    86,     0,    87,     0,    88,     0,
      89,     0,    90,     0,    91,     0,    92,     0,    93,     0,
      94,     0,    95,     0,     0,   155,   156,     0,    69,    96,
       4,     0,    70,    96,     3,     0,    71,    96,     3,     0,
      72,    96,     3,     0,    73,    96,     4,     0,    74,    96,
       3,     0,    75,    96,     3,     0,    76,    96,     3,     0,
      77,    96,     4,     0,    63,     0,    62,     0,     5,     0,
       6,     0,     0,   157,   158,     0,    49,     0,    50,     0,
      51,     0,    52,     0,    53,     0,    54,     0,     7,     0
};

#endif

#if YYDEBUG
/* YYRLINE[YYN] -- source line where rule number YYN was defined. */
static const short yyrline[] =
{
       0,   113,   114,   116,   117,   119,   126,   134,   141,   148,
     155,   161,   167,   173,   179,   182,   190,   196,   202,   206,
     213,   219,   223,   229,   236,   242,   248,   255,   263,   269,
     277,   284,   290,   296,   309,   316,   323,   331,   338,   341,
     347,   355,   361,   367,   374,   375,   377,   378,   383,   388,
     395,   400,   407,   433,   436,   440,   443,   449,   452,   455,
     456,   458,   459,   460,   469,   472,   475,   477,   479,   481,
     483,   484,   486,   489,   492,   495,   496,   498,   500,   506,
     512,   518,   525,   526,   528,   531,   532,   534,   536,   539,
     540,   542,   544,   546,   549,   550,   552,   554,   556,   559,
     560,   562,   567,   569,   575,   576,   578,   580,   582,   584,
     587,   588,   590,   592,   594,   597,   599,   602,   605,   606,
     608,   610,   612,   614,   616,   618,   620,   622,   624,   626,
     628,   630,   632,   634,   636,   638,   640,   642,   644,   646,
     649,   650,   652,   654,   656,   658,   660,   662,   664,   666,
     669,   670,   672,   673,   675,   677,   679,   682,   684,   686,
     688,   690,   692,   694,   696,   698,   700,   702,   704,   706,
     708,   711,   712,   714,   716,   718,   720,   722,   724,   726,
     728,   730,   732,   734,   736,   738,   741,   742,   744,   746,
     748,   750,   752,   754,   757
};
#endif


#if (YYDEBUG) || defined YYERROR_VERBOSE

/* YYTNAME[TOKEN_NUM] -- String name of the token TOKEN_NUM. */
static const char *const yytname[] =
{
  "$", "error", "$undefined.", "REAL", "INTEGER", "QSTRING", "NAME", "NL", 
  "ATTACH", "BUILD", "CALC", "COORD", "CREATE", "DATE", "DOUBLEHELIX", 
  "END", "ENER", "GENERATE", "INIT", "HELIX", "HELP", "LABEL", "MINI", 
  "PRIMARY", "READ", "RECORD", "MODIFY", "SET", "TIME", "WRITE", "STOP", 
  "NSETS", "SAVE", "DIST", "DTAB", "MASS", "VOLUME", "RGYR", "COPY", 
  "ORIENT", "RMS", "STAT", "DERIV", "NFOLD", "OPTIMIZE", "NREPS", "NSTR", 
  "SHIFT", "ONE", "ALL", "BOND", "NBOND", "ANGLE", "TORSION", "ROTBOND", 
  "CARTE", "DIHED", "FRAC", "GRDTOL", "ENETOL", "MAXIT", "MAXFUN", 
  "SILENT", "NOISE", "SYMMETRY", "INTERNAL", "PARAM", "CHARGES", "REF", 
  "PHI", "PHIB", "PHIE", "PHIINC", "PSI", "PSIB", "PSIE", "PSIINC", 
  "RECOPT", "DEBUG", "ID", "TITLE", "RANSEED", "AMBER", "CHARMM", "CFF", 
  "CSD", "DAT", "MM3", "MONO", "PDB", "PF", "PLOT", "SYBYL", "SPF", "PIM", 
  "XYZ", "'='", "'<'", "'>'", "'#'", "'['", "']'", "'('", "':'", "')'", 
  "';'", "'x'", "'y'", "'z'", "'+'", "'-'", "','", "lines", "line", 
  "comstr", "reslist", "residue", "id", "conres", "repeat", "phi", "psi", 
  "ome", "connexion", "conval", "back", "forward", "buildpar", "atoms", 
  "atom", "label", "resno", "coordpars", "coordpar", "enerpars", 
  "enerpar", "genpar", "dhxpars", "dhxpar", "hxpars", "hxpar", "symelems", 
  "symelem", "initsym", "transfs", "transf", "mult", "transl", "sep", 
  "helppar", "minipars", "minipar", "iopars", "iopar", "parformat", 
  "recpars", "recpar", "intpars", "intpar", "term", 0
};
#endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives. */
static const short yyr1[] =
{
       0,   112,   112,   113,   113,   114,   114,   114,   114,   114,
     114,   114,   114,   114,   114,   114,   114,   114,   114,   114,
     114,   114,   114,   114,   114,   114,   114,   114,   114,   114,
     114,   114,   114,   114,   114,   114,   114,   114,   114,   114,
     114,   114,   114,   114,   115,   115,   116,   116,   116,   116,
     116,   116,   116,   117,   118,   119,   120,   121,   122,   123,
     123,   124,   124,   124,   125,   126,   127,   127,   127,   127,
     128,   128,   129,   130,   131,   132,   132,   133,   133,   133,
     133,   133,   134,   134,   135,   136,   136,   136,   136,   137,
     137,   138,   138,   138,   139,   139,   140,   140,   140,   141,
     141,   142,   142,   143,   144,   144,   145,   145,   145,   145,
     146,   146,   146,   146,   146,   147,   147,   148,   149,   149,
     149,   149,   149,   149,   149,   149,   149,   149,   149,   149,
     149,   149,   149,   149,   149,   149,   149,   149,   149,   149,
     150,   150,   151,   151,   151,   151,   151,   151,   151,   151,
     152,   152,   153,   153,   153,   153,   153,   154,   154,   154,
     154,   154,   154,   154,   154,   154,   154,   154,   154,   154,
     154,   155,   155,   156,   156,   156,   156,   156,   156,   156,
     156,   156,   156,   156,   156,   156,   157,   157,   158,   158,
     158,   158,   158,   158,   159
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN. */
static const short yyr2[] =
{
       0,     0,     2,     1,     2,     5,     2,     3,     3,     3,
       2,     2,     2,     2,     2,     3,     1,     2,     1,     2,
       3,     4,     2,     1,     2,     2,     2,     6,     3,     3,
       3,     3,     2,     2,     2,     3,     3,     3,     1,     1,
       3,     3,     2,     1,     0,     2,     1,     4,     3,     1,
       5,     1,     3,     1,     1,     1,     1,     1,     1,     5,
       7,     3,     5,     1,     1,     1,     3,     1,     1,     0,
       0,     2,     2,     1,     1,     0,     2,     1,     2,     1,
       2,     1,     0,     2,     1,     0,     1,     1,     1,     0,
       2,     3,     3,     3,     0,     2,     3,     3,     1,     0,
       2,     7,     1,     1,     0,     2,     2,     2,     2,     1,
       0,     1,     2,     1,     1,     1,     2,     1,     0,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       0,     2,     1,     3,     3,     3,     3,     1,     1,     1,
       0,     2,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     0,     2,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     1,     1,     1,     1,     0,     2,     1,     1,
       1,     1,     1,     1,     1
};

/* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
   doesn't specify something else to do.  Zero means the default is an
   error. */
static const short yydefact[] =
{
       1,     0,    43,   194,     0,    69,     0,    75,     0,    16,
      89,    18,    82,     0,    23,    94,   118,     0,   140,    44,
       0,   171,     0,     0,    39,     0,    38,     2,     0,     3,
      54,     0,     0,    68,    67,     6,    70,    10,    11,    13,
      12,    70,    70,    14,    70,    17,    19,    99,    85,    22,
     119,   120,   121,   122,   123,   124,   125,   126,   127,   128,
     129,   130,   131,   132,   133,   134,   135,   136,   137,   138,
     139,    24,    25,    26,     0,   150,   186,    32,     0,    33,
      70,    34,     0,     0,     0,   150,   186,    42,     4,     0,
       0,     0,     7,     8,     9,    77,    70,    79,    81,    76,
      15,     0,     0,     0,    90,    84,    83,     0,    86,    87,
      88,    20,     0,    98,    95,   149,   142,     0,     0,     0,
       0,   148,   147,   141,    28,    49,    51,     0,    45,    46,
      30,    31,    29,   184,   185,   183,   182,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   172,     0,    35,    37,
      36,    40,    41,    64,     0,    74,     0,    66,    73,    71,
       0,    78,    80,     0,     0,     0,    21,   103,   100,   110,
     102,     0,     0,     0,     0,     0,    53,     0,     0,    55,
       0,   155,   156,   153,   154,   157,   158,   159,   160,   161,
     162,   163,   164,   165,   166,   167,   168,   169,   170,   151,
     152,   188,   189,   190,   191,   192,   193,   187,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     5,
      72,    92,    91,    93,   115,   113,   114,     0,   110,     0,
     109,    96,    97,   143,   144,   145,   146,     0,    48,    52,
     173,   174,   175,   176,   177,   178,   179,   180,   181,     0,
      65,     0,   116,   117,   110,   105,   106,   107,   108,     0,
      47,    27,    59,     0,     0,    50,    56,    63,     0,     0,
     110,     0,    60,     0,    57,    61,   101,     0,    58,    62,
       0,     0
};

static const short yydefgoto[] =
{
       1,    27,    28,    74,   128,   178,    31,   180,   268,   275,
     279,    90,   269,   154,   251,    35,    92,   159,   160,   156,
      43,    99,    46,   106,   111,    45,   104,    49,   114,   107,
     168,   169,   227,   228,   229,   230,   254,    71,    73,   123,
     130,   199,   200,    79,   146,   131,   207,    29
};

static const short yypact[] =
{
  -32768,   133,-32768,-32768,    56,   -17,     4,-32768,   -44,-32768,
  -32768,-32768,-32768,   -43,-32768,-32768,   155,    78,-32768,-32768,
       2,-32768,   -41,   -49,-32768,     5,-32768,-32768,    96,-32768,
  -32768,    14,    38,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
  -32768,-32768,-32768,     6,-32768,   -20,    93,-32768,   -31,    -9,
  -32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
  -32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
  -32768,-32768,-32768,    27,    -2,-32768,-32768,-32768,   131,     3,
  -32768,-32768,   132,   186,   188,-32768,-32768,-32768,-32768,   189,
     190,   191,   192,   192,   192,-32768,-32768,   102,-32768,-32768,
     192,   100,   101,   103,-32768,-32768,-32768,    -6,-32768,-32768,
  -32768,-32768,   104,-32768,-32768,-32768,-32768,   105,   106,   107,
     108,-32768,-32768,-32768,-32768,     1,-32768,   201,-32768,-32768,
      37,   136,-32768,-32768,-32768,-32768,-32768,   110,   111,   112,
     113,   114,   115,   116,   117,   118,-32768,    13,-32768,-32768,
  -32768,    37,   136,-32768,   119,-32768,   210,-32768,-32768,-32768,
     190,   192,-32768,   212,   213,   215,-32768,-32768,-32768,     0,
  -32768,    60,   216,   217,   219,   220,-32768,    56,   123,-32768,
      14,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
  -32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
  -32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,   221,   223,
     224,   225,   226,   228,   229,   230,   231,   138,   232,-32768,
  -32768,-32768,-32768,-32768,   -57,   234,-32768,   127,    -3,   -53,
  -32768,-32768,-32768,-32768,-32768,-32768,-32768,   141,    14,-32768,
  -32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,   237,
  -32768,    10,   -15,-32768,     0,-32768,-32768,-32768,-32768,   143,
  -32768,-32768,-32768,    99,   127,-32768,-32768,-32768,   124,   139,
       8,   239,-32768,   140,-32768,   142,-32768,   242,-32768,-32768,
     246,-32768
};

static const short yypgoto[] =
{
  -32768,-32768,-32768,-32768,-32768,-32768,    71,-32768,-32768,-32768,
  -32768,  -176,-32768,-32768,-32768,-32768,    17,-32768,-32768,    89,
  -32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
  -32768,-32768,  -218,-32768,-32768,-32768,   -14,-32768,-32768,-32768,
     166,-32768,-32768,-32768,-32768,   167,-32768,   -26
};


#define	YYLAST		253


static const short yytable[] =
{
     224,     3,    88,   224,   239,     3,    44,   176,   133,   134,
     255,   224,    80,    75,    32,    33,    85,   108,   109,   158,
      34,    47,    48,   110,   166,   101,   102,   103,   124,    81,
      82,    83,    84,   115,   112,   113,   264,    36,    37,    38,
      39,    40,   181,   182,    95,    96,    97,    98,   129,  -111,
    -111,  -111,   273,   256,   257,   258,    41,    42,    93,    94,
      30,   100,   260,   231,   232,   135,   136,    76,    77,    78,
      86,    87,   137,   138,   139,   140,   141,   142,   143,   144,
     145,   170,   217,   116,    72,   117,   118,   119,   120,   121,
     122,  -112,  -112,  -112,   183,   125,   167,   147,   126,   127,
     177,  -104,   266,     3,   267,   184,   225,   226,  -104,   225,
     226,  -104,  -104,   161,   262,   263,    89,   225,   226,   185,
     186,   187,   188,   189,   190,   191,   192,   193,   194,   195,
     196,   197,   198,   280,    91,   105,   132,   162,   148,     2,
       3,     4,     5,     6,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    50,    51,    52,    53,    54,    55,
      56,    57,    58,    59,    60,    61,    62,    63,    64,    65,
      66,    67,    68,    69,    70,   201,   202,   203,   204,   205,
     206,   149,   150,   153,   155,   157,   163,   164,   158,   165,
     171,   172,   173,   174,   175,   179,   208,   209,   210,   211,
     212,   213,   214,   215,   216,   219,   221,   222,   223,   233,
     234,   238,   218,   235,   236,   240,   241,   242,   243,   271,
     244,   245,   246,   247,   249,   248,   250,   252,   253,   259,
     261,   265,   274,   272,   276,   278,   281,   277,   237,   220,
     270,   151,     0,   152
};

static const short yycheck[] =
{
       3,     7,    28,     3,   180,     7,    50,     6,     5,     6,
     228,     3,    53,    11,    31,    32,    11,    48,    49,     6,
      37,    64,    65,    54,    30,    45,    46,    47,    30,    78,
      79,    80,    81,     6,    43,    44,   254,    33,    34,    35,
      36,    37,     5,     6,    38,    39,    40,    41,    74,   106,
     107,   108,   270,   106,   107,   108,    52,    53,    41,    42,
       4,    44,   238,     3,     4,    62,    63,    65,    66,    67,
      65,    66,    69,    70,    71,    72,    73,    74,    75,    76,
      77,   107,    69,    56,     6,    58,    59,    60,    61,    62,
      63,   106,   107,   108,    57,    97,   102,    80,   100,   101,
      99,   104,     3,     7,     5,    68,   109,   110,   111,   109,
     110,   111,   104,    96,   104,   105,   102,   109,   110,    82,
      83,    84,    85,    86,    87,    88,    89,    90,    91,    92,
      93,    94,    95,     0,    96,    42,     5,    35,     6,     6,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    49,    50,    51,    52,    53,
      54,     5,     4,     4,     4,     4,    96,    96,     6,    96,
      96,    96,    96,    96,    96,     4,    96,    96,    96,    96,
      96,    96,    96,    96,    96,     5,     4,     4,     3,     3,
       3,    98,   103,     4,     4,     4,     3,     3,     3,   105,
       4,     3,     3,     3,    96,     4,     4,     3,   111,    98,
       3,    98,     3,   104,   104,     3,     0,   105,   177,   160,
     264,    85,    -1,    86
};
/* -*-C-*-  Note some compilers choke on comments on `#line' lines.  */
#line 3 "/usr/share/bison/bison.simple"

/* Skeleton output parser for bison,

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002 Free Software
   Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* This is the parser code that is written into each bison parser when
   the %semantic_parser declaration is not specified in the grammar.
   It was written by Richard Stallman by simplifying the hairy parser
   used when %semantic_parser is specified.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

#if ! defined (yyoverflow) || defined (YYERROR_VERBOSE)

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# if YYSTACK_USE_ALLOCA
#  define YYSTACK_ALLOC alloca
# else
#  ifndef YYSTACK_USE_ALLOCA
#   if defined (alloca) || defined (_ALLOCA_H)
#    define YYSTACK_ALLOC alloca
#   else
#    ifdef __GNUC__
#     define YYSTACK_ALLOC __builtin_alloca
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC malloc
#  define YYSTACK_FREE free
# endif
#endif /* ! defined (yyoverflow) || defined (YYERROR_VERBOSE) */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (YYLTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short yyss;
  YYSTYPE yyvs;
# if YYLSP_NEEDED
  YYLTYPE yyls;
# endif
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAX (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# if YYLSP_NEEDED
#  define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE) + sizeof (YYLTYPE))	\
      + 2 * YYSTACK_GAP_MAX)
# else
#  define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE))				\
      + YYSTACK_GAP_MAX)
# endif

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  register YYSIZE_T yyi;		\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (0)
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAX;	\
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (0)

#endif


#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		-2
#define YYEOF		0
#define YYACCEPT	goto yyacceptlab
#define YYABORT 	goto yyabortlab
#define YYERROR		goto yyerrlab1
/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */
#define YYFAIL		goto yyerrlab
#define YYRECOVERING()  (!!yyerrstatus)
#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yychar1 = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { 								\
      yyerror ("syntax error: cannot back up");			\
      YYERROR;							\
    }								\
while (0)

#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Compute the default location (before the actions
   are run).

   When YYLLOC_DEFAULT is run, CURRENT is set the location of the
   first token.  By default, to implement support for ranges, extend
   its range to the last symbol.  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)       	\
   Current.last_line   = Rhs[N].last_line;	\
   Current.last_column = Rhs[N].last_column;
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#if YYPURE
# if YYLSP_NEEDED
#  ifdef YYLEX_PARAM
#   define YYLEX		yylex (&yylval, &yylloc, YYLEX_PARAM)
#  else
#   define YYLEX		yylex (&yylval, &yylloc)
#  endif
# else /* !YYLSP_NEEDED */
#  ifdef YYLEX_PARAM
#   define YYLEX		yylex (&yylval, YYLEX_PARAM)
#  else
#   define YYLEX		yylex (&yylval)
#  endif
# endif /* !YYLSP_NEEDED */
#else /* !YYPURE */
# define YYLEX			yylex ()
#endif /* !YYPURE */


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (0)
/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
#endif /* !YYDEBUG */

/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#if YYMAXDEPTH == 0
# undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif

#ifdef YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  register const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
{
  register char *yyd = yydest;
  register const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif
#endif

#line 315 "/usr/share/bison/bison.simple"


/* The user can define YYPARSE_PARAM as the name of an argument to be passed
   into yyparse.  The argument should have type void *.
   It should actually point to an object.
   Grammar actions can access the variable by casting it
   to the proper pointer type.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
#  define YYPARSE_PARAM_ARG void *YYPARSE_PARAM
#  define YYPARSE_PARAM_DECL
# else
#  define YYPARSE_PARAM_ARG YYPARSE_PARAM
#  define YYPARSE_PARAM_DECL void *YYPARSE_PARAM;
# endif
#else /* !YYPARSE_PARAM */
# define YYPARSE_PARAM_ARG
# define YYPARSE_PARAM_DECL
#endif /* !YYPARSE_PARAM */

/* Prevent warning if -Wstrict-prototypes.  */
#ifdef __GNUC__
# ifdef YYPARSE_PARAM
int yyparse (void *);
# else
int yyparse (void);
# endif
#endif

/* YY_DECL_VARIABLES -- depending whether we use a pure parser,
   variables are global, or local to YYPARSE.  */

#define YY_DECL_NON_LSP_VARIABLES			\
/* The lookahead symbol.  */				\
int yychar;						\
							\
/* The semantic value of the lookahead symbol. */	\
YYSTYPE yylval;						\
							\
/* Number of parse errors so far.  */			\
int yynerrs;

#if YYLSP_NEEDED
# define YY_DECL_VARIABLES			\
YY_DECL_NON_LSP_VARIABLES			\
						\
/* Location data for the lookahead symbol.  */	\
YYLTYPE yylloc;
#else
# define YY_DECL_VARIABLES			\
YY_DECL_NON_LSP_VARIABLES
#endif


/* If nonreentrant, generate the variables here. */

#if !YYPURE
YY_DECL_VARIABLES
#endif  /* !YYPURE */

int
yyparse (YYPARSE_PARAM_ARG)
     YYPARSE_PARAM_DECL
{
  /* If reentrant, generate the variables here. */
#if YYPURE
  YY_DECL_VARIABLES
#endif  /* !YYPURE */

  register int yystate;
  register int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Lookahead token as an internal (translated) token number.  */
  int yychar1 = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack. */
  short	yyssa[YYINITDEPTH];
  short *yyss = yyssa;
  register short *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  register YYSTYPE *yyvsp;

#if YYLSP_NEEDED
  /* The location stack.  */
  YYLTYPE yylsa[YYINITDEPTH];
  YYLTYPE *yyls = yylsa;
  YYLTYPE *yylsp;
#endif

#if YYLSP_NEEDED
# define YYPOPSTACK   (yyvsp--, yyssp--, yylsp--)
#else
# define YYPOPSTACK   (yyvsp--, yyssp--)
#endif

  YYSIZE_T yystacksize = YYINITDEPTH;


  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;
#if YYLSP_NEEDED
  YYLTYPE yyloc;
#endif

  /* When reducing, the number of symbols on the RHS of the reduced
     rule. */
  int yylen;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;
#if YYLSP_NEEDED
  yylsp = yyls;
#endif
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed. so pushing a state here evens the stacks.
     */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyssp >= yyss + yystacksize - 1)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack. Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	short *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  */
# if YYLSP_NEEDED
	YYLTYPE *yyls1 = yyls;
	/* This used to be a conditional around just the two extra args,
	   but that might be undefined if yyoverflow is a macro.  */
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yyls1, yysize * sizeof (*yylsp),
		    &yystacksize);
	yyls = yyls1;
# else
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);
# endif
	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyoverflowlab;
# else
      /* Extend the stack our own way.  */
      if (yystacksize >= YYMAXDEPTH)
	goto yyoverflowlab;
      yystacksize *= 2;
      if (yystacksize > YYMAXDEPTH)
	yystacksize = YYMAXDEPTH;

      {
	short *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyoverflowlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);
# if YYLSP_NEEDED
	YYSTACK_RELOCATE (yyls);
# endif
# undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;
#if YYLSP_NEEDED
      yylsp = yyls + yysize - 1;
#endif

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyssp >= yyss + yystacksize - 1)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* yychar is either YYEMPTY or YYEOF
     or a valid token in external form.  */

  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  /* Convert token to internal form (in yychar1) for indexing tables with */

  if (yychar <= 0)		/* This means end of input. */
    {
      yychar1 = 0;
      yychar = YYEOF;		/* Don't call YYLEX any more */

      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yychar1 = YYTRANSLATE (yychar);

#if YYDEBUG
     /* We have to keep this `#if YYDEBUG', since we use variables
	which are defined only if `YYDEBUG' is set.  */
      if (yydebug)
	{
	  YYFPRINTF (stderr, "Next token is %d (%s",
		     yychar, yytname[yychar1]);
	  /* Give the individual parser a way to print the precise
	     meaning of a token, for further debugging info.  */
# ifdef YYPRINT
	  YYPRINT (stderr, yychar, yylval);
# endif
	  YYFPRINTF (stderr, ")\n");
	}
#endif
    }

  yyn += yychar1;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != yychar1)
    goto yydefault;

  yyn = yytable[yyn];

  /* yyn is what to do for this token type in this state.
     Negative => reduce, -yyn is rule number.
     Positive => shift, yyn is new state.
       New state is final state => don't bother to shift,
       just return success.
     0, or most negative number => error.  */

  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }
  else if (yyn == 0)
    goto yyerrlab;

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */
  YYDPRINTF ((stderr, "Shifting token %d (%s), ",
	      yychar, yytname[yychar1]));

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;
#if YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to the semantic value of
     the lookahead token.  This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];

#if YYLSP_NEEDED
  /* Similarly for the default location.  Let the user run additional
     commands if for instance locations are ranges.  */
  yyloc = yylsp[1-yylen];
  YYLLOC_DEFAULT (yyloc, (yylsp - yylen), yylen);
#endif

#if YYDEBUG
  /* We have to keep this `#if YYDEBUG', since we use variables which
     are defined only if `YYDEBUG' is set.  */
  if (yydebug)
    {
      int yyi;

      YYFPRINTF (stderr, "Reducing via rule %d (line %d), ",
		 yyn, yyrline[yyn]);

      /* Print the symbols being reduced, and their result.  */
      for (yyi = yyprhs[yyn]; yyrhs[yyi] > 0; yyi++)
	YYFPRINTF (stderr, "%s ", yytname[yyrhs[yyi]]);
      YYFPRINTF (stderr, " -> %s\n", yytname[yyr1[yyn]]);
    }
#endif

  switch (yyn) {

case 5:
#line 120 "polypars"
{  putchar('\n');
              putchar('\n');
              macrof = strdup(yylval.string);
              addmacro(polyres, con_b, con_f, no, macrof, phiv, psiv, omev);
              wline('-', MAXLEN);
           ;
    break;}
case 6:
#line 127 "polypars"
{  putchar('\n'); 
              putchar('\n'); 
              tertiary(resno);
              Hydrolysis(da);
              BuildInfo();
              wline('-', MAXLEN);
           ;
    break;}
case 7:
#line 135 "polypars"
{  putchar('\n');
              putchar('\n');
              calc_dist(Aijkl[0], Aijkl[1]);
              inx = 0;
              wline('-', MAXLEN);
           ;
    break;}
case 8:
#line 142 "polypars"
{  putchar('\n');
              putchar('\n');
              calc_angle(Aijkl[0], Aijkl[1], Aijkl[2]);
              inx = 0;
              wline('-', MAXLEN);
           ;
    break;}
case 9:
#line 149 "polypars"
{  putchar('\n');
              putchar('\n');
              calc_tors(Aijkl[0], Aijkl[1], Aijkl[2], Aijkl[3]);
              inx = 0;
              wline('-', MAXLEN);
           ;
    break;}
case 10:
#line 156 "polypars"
{  putchar('\n');
              putchar('\n');
              dist_tabel();
              wline('-', MAXLEN);
           ;
    break;}
case 11:
#line 162 "polypars"
{  putchar('\n');
              putchar('\n');
              calc_mass();
              wline('-', MAXLEN);
           ;
    break;}
case 12:
#line 168 "polypars"
{  putchar('\n');
              putchar('\n');
              calc_rgyr();
              wline('-', MAXLEN);
           ;
    break;}
case 13:
#line 174 "polypars"
{  putchar('\n');
              putchar('\n');
              calc_vol();
              wline('-', MAXLEN);
           ;
    break;}
case 14:
#line 180 "polypars"
{  inx = 0;
           ;
    break;}
case 15:
#line 183 "polypars"
{  putchar('\n');
              putchar('\n');
              add_bond(Aijkl[0],Aijkl[1]);
              calc_dist(Aijkl[0], Aijkl[1]);
              inx = 0;
              wline('-', MAXLEN);
           ;
    break;}
case 16:
#line 191 "polypars"
{  putchar('\n');
              putchar('\n');
              pdate(tid);
              wline('-', MAXLEN);
           ;
    break;}
case 17:
#line 197 "polypars"
{  putchar('\n');
              putchar('\n');
              doubleh(nstr,nreps,shift);
              wline('-', MAXLEN);  
           ;
    break;}
case 18:
#line 203 "polypars"
{  printf("\n\n\tbye - bye\n\n");
              return (0); 
           ;
    break;}
case 19:
#line 207 "polypars"
{  putchar('\n');
              putchar('\n');
              printener(deriv); 
              deriv = 0;
              wline('-', MAXLEN);
           ;
    break;}
case 20:
#line 214 "polypars"
{  putchar('\n');
              putchar('\n');
              Generate(pargenopt); 
              wline('-', MAXLEN);
           ;
    break;}
case 21:
#line 220 "polypars"
{  putchar('\n');
              wline('-', MAXLEN);  
           ;
    break;}
case 22:
#line 224 "polypars"
{  putchar('\n');
              putchar('\n');
              calc_helix(optimize, nfold);
              wline('-', MAXLEN);  
           ;
    break;}
case 23:
#line 230 "polypars"
{  putchar('\n');
              putchar('\n');
              printf("Initializing data structures......\n");
              Initialize();
              wline('-', MAXLEN);
           ;
    break;}
case 24:
#line 237 "polypars"
{  putchar('\n');
              putchar('\n');
              phelp(yylval.cmd);
              wline('-', MAXLEN);
           ;
    break;}
case 25:
#line 243 "polypars"
{  putchar('\n');
              putchar('\n');
              printf("Label: '%s' in line %d\n", yylval.string, lineno);
              wline('-', MAXLEN);
           ;
    break;}
case 26:
#line 249 "polypars"
{  putchar('\n');
              putchar('\n');
              Minimize(parspace,maxfun,maxit,enetol,grdtol,silent);
              parspace = 0;
              wline('-', MAXLEN);
           ;
    break;}
case 27:
#line 256 "polypars"
{  putchar('\n');
              putchar('\n');
              phiv = yylval.real;
              Rotate(Aijkl[0], Aijkl[1], Aijkl[2], Aijkl[3], phiv);
              inx = 0;
              wline('-', MAXLEN);
           ;
    break;}
case 28:
#line 264 "polypars"
{  putchar('\n'); 
              putchar('\n'); 
              res_list(resno);
              wline('-', MAXLEN);
           ;
    break;}
case 29:
#line 270 "polypars"
{  putchar('\n');
              putchar('\n');
              fname = strdup(yylval.string);
              ReadCharges(fname);
              free(fname);
              wline('-', MAXLEN);
           ;
    break;}
case 30:
#line 278 "polypars"
{  putchar('\n');
              putchar('\n');
              ReadCoord(parformat,fname,reference);
              free(fname);
              wline('-', MAXLEN);
           ;
    break;}
case 31:
#line 285 "polypars"
{  putchar('\n');
              putchar('\n');
              printf("Reading internal coordinates......\n");
              wline('-', MAXLEN);
           ;
    break;}
case 32:
#line 291 "polypars"
{  putchar('\n');
              putchar('\n');
              ReadParam();
              wline('-', MAXLEN);
           ;
    break;}
case 33:
#line 297 "polypars"
{  putchar('\n');
              putchar('\n');
              if (psino < 0) 
                 phidriv(phino, phib, phie, phiinc,
                         silent, recopt, fname);
              else
                 phipsidriv(phino, phib, phie, phiinc, 
                            psino, psib, psie, psiinc,
                            silent, recopt, fname);
              free(fname);
              wline('-', MAXLEN);
           ;
    break;}
case 34:
#line 310 "polypars"
{  putchar('\n');
              putchar('\n');
              debug = TRUE;
              printf("Debugging mode on\n");
              wline('-', MAXLEN);
           ;
    break;}
case 35:
#line 317 "polypars"
{  putchar('\n');
              putchar('\n');
              M.id = strdup(yylval.string);
              printf("ID: >%s<\n", M.id);
              wline('-', MAXLEN);
           ;
    break;}
case 36:
#line 324 "polypars"
{  putchar('\n');
              putchar('\n');
              seed = -(yylval.integer);
              printf("The random number generator has been initialized with");
              printf(" a new seed (%d)\n", yylval.integer);
              wline('-', MAXLEN);
           ;
    break;}
case 37:
#line 332 "polypars"
{  putchar('\n');
              putchar('\n');
              M.titel = strdup(yylval.string);
              printf("TITLE: %s\n", M.titel);
              wline('-', MAXLEN);
           ;
    break;}
case 38:
#line 339 "polypars"
{
           ;
    break;}
case 39:
#line 342 "polypars"
{  putchar('\n');
              putchar('\n');
              ptime(tid);
              wline('-', MAXLEN);
           ;
    break;}
case 40:
#line 348 "polypars"
{  putchar('\n');
              putchar('\n');
              WritCoord(parformat, fname, descarte); 
              descarte = TRUE;
              free(fname);
              wline('-', MAXLEN);
           ;
    break;}
case 41:
#line 356 "polypars"
{  putchar('\n');
              putchar('\n');
              list_intcoor(intopt); 
              wline('-', MAXLEN);
           ;
    break;}
case 42:
#line 362 "polypars"
{  putchar('\n');
              putchar('\n');
              printf("\nWriting parameters ........\n"); 
              wline('-', MAXLEN);
           ;
    break;}
case 43:
#line 368 "polypars"
{  putchar('\n');
              putchar('\n');
              yyerror("unknown command:"); 
              wline('-', MAXLEN);
           ;
    break;}
case 47:
#line 379 "polypars"
{ resno++;
             expand(resno,fname,con_b,con_f,thev,phiv,psiv,omev,resno+1,status);
             free(fname);
           ;
    break;}
case 48:
#line 384 "polypars"
{  resno++;
              expand(resno,fname,0,0,999.9,999.9,999.9,999.9,0,status);
              free(fname);
           ;
    break;}
case 49:
#line 389 "polypars"
{  oldstatus=status;
              status++;
              /* if we are in a repeat sequence then save fconnect */
              if ((repno > 0) AND (reptab[repno].fconnect==0))
                 reptab[repno].fconnect = resno;
           ;
    break;}
case 50:
#line 396 "polypars"
{   
              modify_con(resno,polyres,status);
              status=oldstatus;
           ;
    break;}
case 51:
#line 401 "polypars"
{  repno++;
              reptab[repno].fconnect=0;
              reptab[repno].last=0;
              reptab[repno].rep=0;
              reptab[repno].first = resno+1;
           ;
    break;}
case 52:
#line 408 "polypars"
{  
              reptab[repno].last = resno;
              offset = (reptab[repno].last - reptab[repno].first + 1);

              /* first modify the fconnect residue */
              if (reptab[repno].fconnect==0)
                reptab[repno].fconnect=resno;
              modify_con(reptab[repno].fconnect,resno+1,status);

             /* then set the appropiate parameters */
             modify_val(reptab[repno].fconnect,con_b,con_f,thev,phiv,psiv,omev);

              /* loop over number of repeats */
              for (i=1; i<reptab[repno].rep; i++)
                /* loop over offset */
                for (j=1; j<=offset; j++)
                { resno++;
                  copy_res(reptab[repno].first-1+j,resno,i*offset);
                }

       /* if no sidechains to repeat segment connect to the last residue */

              repno--;
           ;
    break;}
case 53:
#line 434 "polypars"
{  fname = strdup(yylval.string); ;
    break;}
case 54:
#line 437 "polypars"
{  polyres = yylval.integer; 
           ;
    break;}
case 55:
#line 441 "polypars"
{  reptab[repno].rep = yylval.integer; ;
    break;}
case 56:
#line 444 "polypars"
{  phiv = 999.9;
              psiv = 999.9;
              omev = 999.9;
              phiv = yylval.real; ;
    break;}
case 57:
#line 450 "polypars"
{  psiv = yylval.real; ;
    break;}
case 58:
#line 453 "polypars"
{  omev = yylval.real; ;
    break;}
case 63:
#line 461 "polypars"
{  phiv = 999.9;
              psiv = 999.9;
              omev = 999.9;
              mapfile = strdup(yylval.string); 
              getconval(mapfile,&phiv,&psiv,&omev);
              free(mapfile);
           ;
    break;}
case 64:
#line 470 "polypars"
{  con_b = yylval.integer; ;
    break;}
case 65:
#line 473 "polypars"
{  con_f = yylval.integer; ;
    break;}
case 66:
#line 476 "polypars"
{  nsets  =  yylval.integer; ;
    break;}
case 67:
#line 478 "polypars"
{ ;
    break;}
case 68:
#line 480 "polypars"
{ ;
    break;}
case 72:
#line 487 "polypars"
{  Aijkl[inx++] = getAno(no, lab); ;
    break;}
case 73:
#line 490 "polypars"
{  strcpy(lab, yyval.string); ;
    break;}
case 74:
#line 493 "polypars"
{  no = yyval.integer; ;
    break;}
case 77:
#line 499 "polypars"
{  Coord_copy(); ;
    break;}
case 78:
#line 501 "polypars"
{  putchar('\n');
              putchar('\n');
              Coord_orient(Aijkl[0],Aijkl[1],Aijkl[2]);
              wline('-', MAXLEN);
           ;
    break;}
case 79:
#line 507 "polypars"
{  putchar('\n');
              putchar('\n');
              Coord_rms(FALSE);
              wline('-', MAXLEN);
           ;
    break;}
case 80:
#line 513 "polypars"
{  putchar('\n');
              putchar('\n');
              Coord_rms(TRUE);
              wline('-', MAXLEN);
           ;
    break;}
case 81:
#line 519 "polypars"
{  putchar('\n');
              putchar('\n');
              Coord_stat(); 
              wline('-', MAXLEN);
           ;
    break;}
case 84:
#line 529 "polypars"
{  deriv = 1; ;
    break;}
case 85:
#line 531 "polypars"
{  pargenopt = 1; ;
    break;}
case 86:
#line 533 "polypars"
{  pargenopt = 0; ;
    break;}
case 87:
#line 535 "polypars"
{  pargenopt = 1; ;
    break;}
case 88:
#line 537 "polypars"
{  pargenopt = 11; ;
    break;}
case 91:
#line 543 "polypars"
{  nstr = yylval.integer; ;
    break;}
case 92:
#line 545 "polypars"
{  nreps = yylval.integer; ;
    break;}
case 93:
#line 547 "polypars"
{  shift = yylval.real; ;
    break;}
case 96:
#line 553 "polypars"
{  nfold = yylval.real; ;
    break;}
case 97:
#line 555 "polypars"
{  nfold = (double) yylval.integer; ;
    break;}
case 98:
#line 557 "polypars"
{  optimize = TRUE; ;
    break;}
case 101:
#line 563 "polypars"
{  if (rnat == 0)
                 rnat = M.nat;
              symneighbor(rnat, TM); 
           ;
    break;}
case 103:
#line 570 "polypars"
{  initmat(TM, 0.0); 
              TM[3][3] = 1.0;
              symno = 0;
           ;
    break;}
case 106:
#line 579 "polypars"
{  TM[symno][0] = xmult; ;
    break;}
case 107:
#line 581 "polypars"
{  TM[symno][1] = xmult; ;
    break;}
case 108:
#line 583 "polypars"
{  TM[symno][2] = xmult; ;
    break;}
case 109:
#line 585 "polypars"
{  TM[symno][3] = xtransl; ;
    break;}
case 110:
#line 587 "polypars"
{  xmult = 1.0; ;
    break;}
case 111:
#line 589 "polypars"
{  xmult = yylval.real; ;
    break;}
case 112:
#line 591 "polypars"
{  xmult = yylval.real; ;
    break;}
case 113:
#line 593 "polypars"
{  xmult = 1.0; ;
    break;}
case 114:
#line 595 "polypars"
{  xmult = -1.0; ;
    break;}
case 115:
#line 598 "polypars"
{  xtransl = yylval.real; ;
    break;}
case 116:
#line 600 "polypars"
{  xtransl = yylval.real; ;
    break;}
case 117:
#line 603 "polypars"
{  symno++; ;
    break;}
case 119:
#line 607 "polypars"
{ ;
    break;}
case 120:
#line 609 "polypars"
{ ;
    break;}
case 121:
#line 611 "polypars"
{ ;
    break;}
case 122:
#line 613 "polypars"
{ ;
    break;}
case 123:
#line 615 "polypars"
{ ;
    break;}
case 124:
#line 617 "polypars"
{ ;
    break;}
case 125:
#line 619 "polypars"
{ ;
    break;}
case 126:
#line 621 "polypars"
{ ;
    break;}
case 127:
#line 623 "polypars"
{ ;
    break;}
case 128:
#line 625 "polypars"
{ ;
    break;}
case 129:
#line 627 "polypars"
{ ;
    break;}
case 130:
#line 629 "polypars"
{ ;
    break;}
case 131:
#line 631 "polypars"
{ ;
    break;}
case 132:
#line 633 "polypars"
{ ;
    break;}
case 133:
#line 635 "polypars"
{ ;
    break;}
case 134:
#line 637 "polypars"
{ ;
    break;}
case 135:
#line 639 "polypars"
{ ;
    break;}
case 136:
#line 641 "polypars"
{ ;
    break;}
case 137:
#line 643 "polypars"
{ ;
    break;}
case 138:
#line 645 "polypars"
{ ;
    break;}
case 139:
#line 647 "polypars"
{ ;
    break;}
case 142:
#line 653 "polypars"
{  parspace = 1; ;
    break;}
case 143:
#line 655 "polypars"
{  grdtol =  yylval.real; ;
    break;}
case 144:
#line 657 "polypars"
{  enetol =  yylval.real; ;
    break;}
case 145:
#line 659 "polypars"
{  maxit  =  yylval.integer; ;
    break;}
case 146:
#line 661 "polypars"
{  maxfun =  yylval.integer; ;
    break;}
case 147:
#line 663 "polypars"
{  silent = FALSE; ;
    break;}
case 148:
#line 665 "polypars"
{  silent = TRUE; ;
    break;}
case 149:
#line 667 "polypars"
{ yyerror("unknown subcommand: "); ;
    break;}
case 153:
#line 674 "polypars"
{  descarte = FALSE; ;
    break;}
case 154:
#line 676 "polypars"
{  reference = TRUE; ;
    break;}
case 155:
#line 678 "polypars"
{  fname = strdup(yylval.string); ;
    break;}
case 156:
#line 680 "polypars"
{ yyerror("unknown subcommand: "); ;
    break;}
case 157:
#line 683 "polypars"
{ parformat = 0; ;
    break;}
case 158:
#line 685 "polypars"
{ parformat = 1; ;
    break;}
case 159:
#line 687 "polypars"
{ parformat = 2; ;
    break;}
case 160:
#line 689 "polypars"
{ parformat = 3; ;
    break;}
case 161:
#line 691 "polypars"
{ parformat = 4; ;
    break;}
case 162:
#line 693 "polypars"
{ parformat = 5; ;
    break;}
case 163:
#line 695 "polypars"
{ parformat = 6; ;
    break;}
case 164:
#line 697 "polypars"
{ parformat = 7; ;
    break;}
case 165:
#line 699 "polypars"
{ parformat = 8; ;
    break;}
case 166:
#line 701 "polypars"
{ parformat = 9; ;
    break;}
case 167:
#line 703 "polypars"
{ parformat = 10; ;
    break;}
case 168:
#line 705 "polypars"
{ parformat = 11; ;
    break;}
case 169:
#line 707 "polypars"
{ parformat = 12; ;
    break;}
case 170:
#line 709 "polypars"
{ parformat = 13; ;
    break;}
case 173:
#line 715 "polypars"
{  phino = yylval.integer; ;
    break;}
case 174:
#line 717 "polypars"
{  phib  = yylval.real; ;
    break;}
case 175:
#line 719 "polypars"
{  phie  = yylval.real; ;
    break;}
case 176:
#line 721 "polypars"
{  phiinc= yylval.real; ;
    break;}
case 177:
#line 723 "polypars"
{  psino = yylval.integer; ;
    break;}
case 178:
#line 725 "polypars"
{  psib  = yylval.real; ;
    break;}
case 179:
#line 727 "polypars"
{  psie  = yylval.real; ;
    break;}
case 180:
#line 729 "polypars"
{  psiinc= yylval.real; ;
    break;}
case 181:
#line 731 "polypars"
{  recopt = yylval.integer; ;
    break;}
case 182:
#line 733 "polypars"
{  silent = FALSE; ;
    break;}
case 183:
#line 735 "polypars"
{  silent = TRUE; ;
    break;}
case 184:
#line 737 "polypars"
{  fname = strdup(yylval.string); ;
    break;}
case 185:
#line 739 "polypars"
{ yyerror("unknown subcommand: "); ;
    break;}
case 186:
#line 741 "polypars"
{ intopt = 1; ;
    break;}
case 188:
#line 745 "polypars"
{ intopt = 2; ;
    break;}
case 189:
#line 747 "polypars"
{ intopt = 1; ;
    break;}
case 190:
#line 749 "polypars"
{ intopt = 1; ;
    break;}
case 191:
#line 751 "polypars"
{ intopt = 1; ;
    break;}
case 192:
#line 753 "polypars"
{ intopt = 1; ;
    break;}
case 193:
#line 755 "polypars"
{ intopt = 0; ;
    break;}
case 194:
#line 758 "polypars"
{  lineno++;
              if (lineno > 10000) lineno=1;
              printf("\nPOLYS:%d > ", lineno);
           ;
    break;}
}

#line 705 "/usr/share/bison/bison.simple"


  yyvsp -= yylen;
  yyssp -= yylen;
#if YYLSP_NEEDED
  yylsp -= yylen;
#endif

#if YYDEBUG
  if (yydebug)
    {
      short *yyssp1 = yyss - 1;
      YYFPRINTF (stderr, "state stack now");
      while (yyssp1 != yyssp)
	YYFPRINTF (stderr, " %d", *++yyssp1);
      YYFPRINTF (stderr, "\n");
    }
#endif

  *++yyvsp = yyval;
#if YYLSP_NEEDED
  *++yylsp = yyloc;
#endif

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTBASE] + *yyssp;
  if (yystate >= 0 && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTBASE];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;

#ifdef YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (yyn > YYFLAG && yyn < YYLAST)
	{
	  YYSIZE_T yysize = 0;
	  char *yymsg;
	  int yyx, yycount;

	  yycount = 0;
	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  for (yyx = yyn < 0 ? -yyn : 0;
	       yyx < (int) (sizeof (yytname) / sizeof (char *)); yyx++)
	    if (yycheck[yyx + yyn] == yyx)
	      yysize += yystrlen (yytname[yyx]) + 15, yycount++;
	  yysize += yystrlen ("parse error, unexpected ") + 1;
	  yysize += yystrlen (yytname[YYTRANSLATE (yychar)]);
	  yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg != 0)
	    {
	      char *yyp = yystpcpy (yymsg, "parse error, unexpected ");
	      yyp = yystpcpy (yyp, yytname[YYTRANSLATE (yychar)]);

	      if (yycount < 5)
		{
		  yycount = 0;
		  for (yyx = yyn < 0 ? -yyn : 0;
		       yyx < (int) (sizeof (yytname) / sizeof (char *));
		       yyx++)
		    if (yycheck[yyx + yyn] == yyx)
		      {
			const char *yyq = ! yycount ? ", expecting " : " or ";
			yyp = yystpcpy (yyp, yyq);
			yyp = yystpcpy (yyp, yytname[yyx]);
			yycount++;
		      }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    yyerror ("parse error; also virtual memory exhausted");
	}
      else
#endif /* defined (YYERROR_VERBOSE) */
	yyerror ("parse error");
    }
  goto yyerrlab1;


/*--------------------------------------------------.
| yyerrlab1 -- error raised explicitly by an action |
`--------------------------------------------------*/
yyerrlab1:
  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      /* return failure if at end of input */
      if (yychar == YYEOF)
	YYABORT;
      YYDPRINTF ((stderr, "Discarding token %d (%s).\n",
		  yychar, yytname[yychar1]));
      yychar = YYEMPTY;
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */

  yyerrstatus = 3;		/* Each real token shifted decrements this */

  goto yyerrhandle;


/*-------------------------------------------------------------------.
| yyerrdefault -- current state does not do anything special for the |
| error token.                                                       |
`-------------------------------------------------------------------*/
yyerrdefault:
#if 0
  /* This is wrong; only states that explicitly want error tokens
     should shift them.  */

  /* If its default is to accept any token, ok.  Otherwise pop it.  */
  yyn = yydefact[yystate];
  if (yyn)
    goto yydefault;
#endif


/*---------------------------------------------------------------.
| yyerrpop -- pop the current state because it cannot handle the |
| error token                                                    |
`---------------------------------------------------------------*/
yyerrpop:
  if (yyssp == yyss)
    YYABORT;
  yyvsp--;
  yystate = *--yyssp;
#if YYLSP_NEEDED
  yylsp--;
#endif

#if YYDEBUG
  if (yydebug)
    {
      short *yyssp1 = yyss - 1;
      YYFPRINTF (stderr, "Error: state stack now");
      while (yyssp1 != yyssp)
	YYFPRINTF (stderr, " %d", *++yyssp1);
      YYFPRINTF (stderr, "\n");
    }
#endif

/*--------------.
| yyerrhandle.  |
`--------------*/
yyerrhandle:
  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yyerrdefault;

  yyn += YYTERROR;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != YYTERROR)
    goto yyerrdefault;

  yyn = yytable[yyn];
  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	goto yyerrpop;
      yyn = -yyn;
      goto yyreduce;
    }
  else if (yyn == 0)
    goto yyerrpop;

  if (yyn == YYFINAL)
    YYACCEPT;

  YYDPRINTF ((stderr, "Shifting error token, "));

  *++yyvsp = yylval;
#if YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

/*---------------------------------------------.
| yyoverflowab -- parser overflow comes here.  |
`---------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}
#line 763 "polypars"


