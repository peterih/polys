
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.4.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Copy the first part of user declarations.  */

/* Line 189 of yacc.c  */
#line 9 "polypars"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "polys.h"

/* parser variables needed for the subcommands */
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

/* parser variables needed for the BUILD command */
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

/* External variabels  needed for the parser */
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


/* Line 189 of yacc.c  */
#line 148 "polypars.tab.c"

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     REAL = 258,
     INTEGER = 259,
     QSTRING = 260,
     NAME = 261,
     NL = 262,
     ATTACH = 263,
     BUILD = 264,
     CALC = 265,
     COORD = 266,
     CREATE = 267,
     DATE = 268,
     DOUBLEHELIX = 269,
     END = 270,
     ENER = 271,
     GENERATE = 272,
     INIT = 273,
     HELIX = 274,
     HELP = 275,
     LABEL = 276,
     MINI = 277,
     PRIMARY = 278,
     READ = 279,
     RECORD = 280,
     MODIFY = 281,
     SET = 282,
     TIME = 283,
     WRITE = 284,
     STOP = 285,
     NSETS = 286,
     SAVE = 287,
     DIST = 288,
     DTAB = 289,
     MASS = 290,
     VOLUME = 291,
     RGYR = 292,
     COPY = 293,
     ORIENT = 294,
     RMS = 295,
     STAT = 296,
     DERIV = 297,
     NFOLD = 298,
     OPTIMIZE = 299,
     NREPS = 300,
     NSTR = 301,
     SHIFT = 302,
     ONE = 303,
     ALL = 304,
     BOND = 305,
     NBOND = 306,
     ANGLE = 307,
     TORSION = 308,
     ROTBOND = 309,
     CARTE = 310,
     DIHED = 311,
     FRAC = 312,
     GRDTOL = 313,
     ENETOL = 314,
     MAXIT = 315,
     MAXFUN = 316,
     SILENT = 317,
     NOISE = 318,
     SYMMETRY = 319,
     INTERNAL = 320,
     PARAM = 321,
     CHARGES = 322,
     REF = 323,
     PHI = 324,
     PHIB = 325,
     PHIE = 326,
     PHIINC = 327,
     PSI = 328,
     PSIB = 329,
     PSIE = 330,
     PSIINC = 331,
     RECOPT = 332,
     DEBUG = 333,
     ID = 334,
     TITLE = 335,
     RANSEED = 336,
     AMBER = 337,
     CHARMM = 338,
     CFF = 339,
     CSD = 340,
     DAT = 341,
     MM3 = 342,
     MONO = 343,
     PDB = 344,
     PF = 345,
     PLOT = 346,
     SYBYL = 347,
     SPF = 348,
     PIM = 349,
     XYZ = 350
   };
#endif



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 82 "polypars"
 double   real;
         int      integer;
         char     *string;
         int      cmd;
       


/* Line 214 of yacc.c  */
#line 287 "polypars.tab.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 299 "polypars.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   241

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  112
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  49
/* YYNRULES -- Number of rules.  */
#define YYNRULES  194
/* YYNRULES -- Number of states.  */
#define YYNSTATES  281

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   350

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
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
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     4,     7,     9,    12,    18,    21,    25,
      29,    33,    36,    39,    42,    45,    48,    52,    54,    57,
      59,    62,    66,    71,    74,    76,    79,    82,    85,    92,
      96,   100,   104,   108,   111,   114,   117,   121,   125,   129,
     131,   133,   137,   141,   144,   146,   147,   150,   152,   157,
     161,   163,   169,   171,   175,   177,   179,   181,   183,   185,
     187,   193,   201,   205,   211,   213,   215,   217,   221,   223,
     225,   226,   227,   230,   233,   235,   237,   238,   241,   243,
     246,   248,   251,   253,   254,   257,   259,   260,   262,   264,
     266,   267,   270,   274,   278,   282,   283,   286,   290,   294,
     296,   297,   300,   308,   310,   312,   313,   316,   319,   322,
     325,   327,   329,   332,   334,   336,   338,   341,   343,   344,
     346,   348,   350,   352,   354,   356,   358,   360,   362,   364,
     366,   368,   370,   372,   374,   376,   378,   380,   382,   384,
     386,   387,   390,   392,   396,   400,   404,   408,   410,   412,
     414,   415,   418,   420,   422,   424,   426,   428,   430,   432,
     434,   436,   438,   440,   442,   444,   446,   448,   450,   452,
     454,   456,   457,   460,   464,   468,   472,   476,   480,   484,
     488,   492,   496,   498,   500,   502,   504,   505,   508,   510,
     512,   514,   516,   518,   520
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
     113,     0,    -1,    -1,   113,   114,    -1,   160,    -1,   115,
     160,    -1,     8,   119,   124,   132,     5,    -1,     9,   128,
      -1,    10,    33,   129,    -1,    10,    52,   129,    -1,    10,
      53,   129,    -1,    10,    34,    -1,    10,    35,    -1,    10,
      37,    -1,    10,    36,    -1,    11,   133,    -1,    12,    50,
     129,    -1,    13,    -1,    14,   138,    -1,    15,    -1,    16,
     135,    -1,    17,    65,   137,    -1,    17,    64,   142,    30,
      -1,    19,   140,    -1,    18,    -1,    20,   150,    -1,    21,
       6,    -1,    22,   151,    -1,    26,    53,   129,    69,    96,
       3,    -1,    23,   116,    30,    -1,    24,    67,     5,    -1,
      24,    11,   153,    -1,    24,    65,   158,    -1,    24,    66,
      -1,    25,   156,    -1,    27,    78,    -1,    27,    79,     6,
      -1,    27,    81,     4,    -1,    27,    80,     5,    -1,    30,
      -1,    28,    -1,    29,    11,   153,    -1,    29,    65,   158,
      -1,    29,    66,    -1,     6,    -1,    -1,   116,   117,    -1,
     160,    -1,    97,   118,    98,   124,    -1,    97,   118,    98,
      -1,    97,    -1,    97,    99,   119,    98,    98,    -1,   100,
      -1,   101,   120,   124,    -1,     6,    -1,     4,    -1,     4,
      -1,     3,    -1,     3,    -1,     3,    -1,   102,   126,   103,
     127,   104,    -1,   102,   126,   103,   127,   105,   125,   104,
      -1,   121,   105,   122,    -1,   121,   105,   122,   105,   123,
      -1,     5,    -1,     4,    -1,     4,    -1,    31,    96,     4,
      -1,    37,    -1,    32,    -1,    -1,    -1,   129,   130,    -1,
     131,   132,    -1,     6,    -1,     4,    -1,    -1,   133,   134,
      -1,    38,    -1,    39,   129,    -1,    40,    -1,    40,    35,
      -1,    41,    -1,    -1,   135,   136,    -1,    42,    -1,    -1,
      48,    -1,    49,    -1,    54,    -1,    -1,   138,   139,    -1,
      46,    96,     4,    -1,    45,    96,     4,    -1,    47,    96,
       3,    -1,    -1,   140,   141,    -1,    43,    96,     3,    -1,
      43,    96,     4,    -1,    44,    -1,    -1,   142,   143,    -1,
     144,   145,   149,   145,   149,   145,   104,    -1,   160,    -1,
     102,    -1,    -1,   146,   145,    -1,   147,   106,    -1,   147,
     107,    -1,   147,   108,    -1,   148,    -1,     3,    -1,   109,
       3,    -1,   109,    -1,   110,    -1,     3,    -1,   109,     3,
      -1,   111,    -1,    -1,     9,    -1,    10,    -1,    11,    -1,
      12,    -1,    13,    -1,    14,    -1,    15,    -1,    16,    -1,
      17,    -1,    18,    -1,    19,    -1,    20,    -1,    21,    -1,
      22,    -1,    23,    -1,    24,    -1,    25,    -1,    26,    -1,
      27,    -1,    28,    -1,    29,    -1,    -1,   151,   152,    -1,
      56,    -1,    58,    96,     3,    -1,    59,    96,     3,    -1,
      60,    96,     4,    -1,    61,    96,     4,    -1,    63,    -1,
      62,    -1,     6,    -1,    -1,   153,   154,    -1,   155,    -1,
      57,    -1,    68,    -1,     5,    -1,     6,    -1,    82,    -1,
      83,    -1,    84,    -1,    85,    -1,    86,    -1,    87,    -1,
      88,    -1,    89,    -1,    90,    -1,    91,    -1,    92,    -1,
      93,    -1,    94,    -1,    95,    -1,    -1,   156,   157,    -1,
      69,    96,     4,    -1,    70,    96,     3,    -1,    71,    96,
       3,    -1,    72,    96,     3,    -1,    73,    96,     4,    -1,
      74,    96,     3,    -1,    75,    96,     3,    -1,    76,    96,
       3,    -1,    77,    96,     4,    -1,    63,    -1,    62,    -1,
       5,    -1,     6,    -1,    -1,   158,   159,    -1,    49,    -1,
      50,    -1,    51,    -1,    52,    -1,    53,    -1,    54,    -1,
       7,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   107,   107,   108,   110,   111,   113,   120,   128,   135,
     142,   149,   155,   161,   167,   173,   176,   184,   190,   196,
     200,   207,   213,   217,   223,   230,   236,   242,   249,   257,
     263,   271,   278,   284,   290,   303,   310,   317,   325,   332,
     335,   341,   349,   355,   361,   368,   369,   371,   372,   377,
     382,   389,   394,   401,   427,   430,   434,   437,   443,   446,
     449,   450,   452,   453,   454,   463,   466,   469,   471,   473,
     475,   477,   478,   480,   483,   486,   489,   490,   492,   494,
     500,   506,   512,   519,   520,   522,   525,   526,   528,   530,
     533,   534,   536,   538,   540,   543,   544,   546,   548,   550,
     553,   554,   556,   561,   563,   569,   570,   572,   574,   576,
     578,   581,   583,   585,   587,   590,   592,   595,   598,   599,
     601,   603,   605,   607,   609,   611,   613,   615,   617,   619,
     621,   623,   625,   627,   629,   631,   633,   635,   637,   639,
     642,   643,   645,   647,   649,   651,   653,   655,   657,   659,
     662,   663,   665,   666,   668,   670,   672,   675,   677,   679,
     681,   683,   685,   687,   689,   691,   693,   695,   697,   699,
     701,   704,   705,   707,   709,   711,   713,   715,   717,   719,
     721,   723,   725,   727,   729,   731,   734,   735,   737,   739,
     741,   743,   745,   747,   750
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "REAL", "INTEGER", "QSTRING", "NAME",
  "NL", "ATTACH", "BUILD", "CALC", "COORD", "CREATE", "DATE",
  "DOUBLEHELIX", "END", "ENER", "GENERATE", "INIT", "HELIX", "HELP",
  "LABEL", "MINI", "PRIMARY", "READ", "RECORD", "MODIFY", "SET", "TIME",
  "WRITE", "STOP", "NSETS", "SAVE", "DIST", "DTAB", "MASS", "VOLUME",
  "RGYR", "COPY", "ORIENT", "RMS", "STAT", "DERIV", "NFOLD", "OPTIMIZE",
  "NREPS", "NSTR", "SHIFT", "ONE", "ALL", "BOND", "NBOND", "ANGLE",
  "TORSION", "ROTBOND", "CARTE", "DIHED", "FRAC", "GRDTOL", "ENETOL",
  "MAXIT", "MAXFUN", "SILENT", "NOISE", "SYMMETRY", "INTERNAL", "PARAM",
  "CHARGES", "REF", "PHI", "PHIB", "PHIE", "PHIINC", "PSI", "PSIB", "PSIE",
  "PSIINC", "RECOPT", "DEBUG", "ID", "TITLE", "RANSEED", "AMBER", "CHARMM",
  "CFF", "CSD", "DAT", "MM3", "MONO", "PDB", "PF", "PLOT", "SYBYL", "SPF",
  "PIM", "XYZ", "'='", "'<'", "'>'", "'#'", "'['", "']'", "'('", "':'",
  "')'", "';'", "'x'", "'y'", "'z'", "'+'", "'-'", "','", "$accept",
  "lines", "line", "comstr", "reslist", "residue", "id", "conres",
  "repeat", "phi", "psi", "ome", "connexion", "conval", "back", "forward",
  "buildpar", "atoms", "atom", "label", "resno", "coordpars", "coordpar",
  "enerpars", "enerpar", "genpar", "dhxpars", "dhxpar", "hxpars", "hxpar",
  "symelems", "symelem", "initsym", "transfs", "transf", "mult", "transl",
  "sep", "helppar", "minipars", "minipar", "iopars", "iopar", "parformat",
  "recpars", "recpar", "intpars", "intpar", "term", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,    61,    60,    62,    35,
      91,    93,    40,    58,    41,    59,   120,   121,   122,    43,
      45,    44
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,   112,   113,   113,   114,   114,   115,   115,   115,   115,
     115,   115,   115,   115,   115,   115,   115,   115,   115,   115,
     115,   115,   115,   115,   115,   115,   115,   115,   115,   115,
     115,   115,   115,   115,   115,   115,   115,   115,   115,   115,
     115,   115,   115,   115,   115,   116,   116,   117,   117,   117,
     117,   117,   117,   117,   118,   119,   120,   121,   122,   123,
     124,   124,   125,   125,   125,   126,   127,   128,   128,   128,
     128,   129,   129,   130,   131,   132,   133,   133,   134,   134,
     134,   134,   134,   135,   135,   136,   137,   137,   137,   137,
     138,   138,   139,   139,   139,   140,   140,   141,   141,   141,
     142,   142,   143,   143,   144,   145,   145,   146,   146,   146,
     146,   147,   147,   147,   147,   148,   148,   149,   150,   150,
     150,   150,   150,   150,   150,   150,   150,   150,   150,   150,
     150,   150,   150,   150,   150,   150,   150,   150,   150,   150,
     151,   151,   152,   152,   152,   152,   152,   152,   152,   152,
     153,   153,   154,   154,   154,   154,   154,   155,   155,   155,
     155,   155,   155,   155,   155,   155,   155,   155,   155,   155,
     155,   156,   156,   157,   157,   157,   157,   157,   157,   157,
     157,   157,   157,   157,   157,   157,   158,   158,   159,   159,
     159,   159,   159,   159,   160
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     1,     2,     5,     2,     3,     3,
       3,     2,     2,     2,     2,     2,     3,     1,     2,     1,
       2,     3,     4,     2,     1,     2,     2,     2,     6,     3,
       3,     3,     3,     2,     2,     2,     3,     3,     3,     1,
       1,     3,     3,     2,     1,     0,     2,     1,     4,     3,
       1,     5,     1,     3,     1,     1,     1,     1,     1,     1,
       5,     7,     3,     5,     1,     1,     1,     3,     1,     1,
       0,     0,     2,     2,     1,     1,     0,     2,     1,     2,
       1,     2,     1,     0,     2,     1,     0,     1,     1,     1,
       0,     2,     3,     3,     3,     0,     2,     3,     3,     1,
       0,     2,     7,     1,     1,     0,     2,     2,     2,     2,
       1,     1,     2,     1,     1,     1,     2,     1,     0,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       0,     2,     1,     3,     3,     3,     3,     1,     1,     1,
       0,     2,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     0,     2,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     1,     1,     1,     1,     0,     2,     1,     1,
       1,     1,     1,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,    44,   194,     0,    70,     0,    76,     0,
      17,    90,    19,    83,     0,    24,    95,   118,     0,   140,
      45,     0,   171,     0,     0,    40,     0,    39,     3,     0,
       4,    55,     0,     0,    69,    68,     7,    71,    11,    12,
      14,    13,    71,    71,    15,    71,    18,    20,   100,    86,
      23,   119,   120,   121,   122,   123,   124,   125,   126,   127,
     128,   129,   130,   131,   132,   133,   134,   135,   136,   137,
     138,   139,    25,    26,    27,     0,   150,   186,    33,     0,
      34,    71,    35,     0,     0,     0,   150,   186,    43,     5,
       0,     0,     0,     8,     9,    10,    78,    71,    80,    82,
      77,    16,     0,     0,     0,    91,    85,    84,     0,    87,
      88,    89,    21,     0,    99,    96,   149,   142,     0,     0,
       0,     0,   148,   147,   141,    29,    50,    52,     0,    46,
      47,    31,    32,    30,   184,   185,   183,   182,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   172,     0,    36,
      38,    37,    41,    42,    65,     0,    75,     0,    67,    74,
      72,     0,    79,    81,     0,     0,     0,    22,   104,   101,
     105,   103,     0,     0,     0,     0,     0,    54,     0,     0,
      56,     0,   155,   156,   153,   154,   157,   158,   159,   160,
     161,   162,   163,   164,   165,   166,   167,   168,   169,   170,
     151,   152,   188,   189,   190,   191,   192,   193,   187,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       6,    73,    93,    92,    94,   115,   113,   114,     0,   105,
       0,   110,    97,    98,   143,   144,   145,   146,     0,    49,
      53,   173,   174,   175,   176,   177,   178,   179,   180,   181,
       0,    66,     0,   116,   117,   105,   106,   107,   108,   109,
       0,    48,    28,    60,     0,     0,    51,    57,    64,     0,
       0,   105,     0,    61,     0,    58,    62,   102,     0,    59,
      63
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     1,    28,    29,    75,   129,   179,    32,   181,   269,
     276,   280,    91,   270,   155,   252,    36,    93,   160,   161,
     157,    44,   100,    47,   107,   112,    46,   105,    50,   115,
     108,   169,   170,   228,   229,   230,   231,   255,    72,    74,
     124,   131,   200,   201,    80,   147,   132,   208,    30
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -220
static const yytype_int16 yypact[] =
{
    -220,   122,  -220,  -220,  -220,     7,     8,   -15,  -220,   -36,
    -220,  -220,  -220,  -220,    40,  -220,  -220,   144,    11,  -220,
    -220,    -7,  -220,     4,   -37,  -220,     5,  -220,  -220,    16,
    -220,  -220,   -69,   -50,  -220,  -220,  -220,  -220,  -220,  -220,
    -220,  -220,  -220,  -220,    15,  -220,   -20,    43,  -220,   -19,
      49,  -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,
    -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,
    -220,  -220,  -220,  -220,    28,    -2,  -220,  -220,  -220,    45,
       3,  -220,  -220,    62,    77,    93,  -220,  -220,  -220,  -220,
     178,   179,   180,   181,   181,   181,  -220,  -220,   150,  -220,
    -220,   181,    90,    92,    94,  -220,  -220,  -220,    -6,  -220,
    -220,  -220,  -220,    95,  -220,  -220,  -220,  -220,    96,    97,
      98,    99,  -220,  -220,  -220,  -220,     1,  -220,   185,  -220,
    -220,    26,   125,  -220,  -220,  -220,  -220,  -220,   100,   101,
     102,   103,   104,   105,   106,   107,   108,  -220,     0,  -220,
    -220,  -220,    26,   125,  -220,   109,  -220,   200,  -220,  -220,
    -220,   179,   181,  -220,   202,   203,   205,  -220,  -220,  -220,
      -3,  -220,   123,   206,   207,   209,   210,  -220,     7,   113,
    -220,   -69,  -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,
    -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,
    -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,   211,
     213,   214,   215,   216,   218,   219,   220,   221,   128,   222,
    -220,  -220,  -220,  -220,  -220,   -44,   224,  -220,   117,    -3,
      -5,  -220,  -220,  -220,  -220,  -220,  -220,  -220,   121,   -69,
    -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,
     226,  -220,    76,    17,  -220,    -3,  -220,  -220,  -220,  -220,
     132,  -220,  -220,  -220,    44,   117,  -220,  -220,  -220,   126,
     129,    -3,   229,  -220,   130,  -220,   131,  -220,   232,  -220,
    -220
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -220,  -220,  -220,  -220,  -220,  -220,  -220,    59,  -220,  -220,
    -220,  -220,  -178,  -220,  -220,  -220,  -220,   -30,  -220,  -220,
      78,  -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,
    -220,  -220,  -220,  -219,  -220,  -220,  -220,   -25,  -220,  -220,
    -220,   152,  -220,  -220,  -220,  -220,   154,  -220,   -27
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -113
static const yytype_int16 yytable[] =
{
     225,     4,    89,   240,    76,     4,   159,   177,   134,   135,
     256,    31,    94,    95,    45,   101,    86,    73,    37,    38,
      39,    40,    41,     4,   167,   102,   103,   104,   125,   109,
     110,   182,   183,    90,   116,   111,   265,    42,    43,    33,
      34,    82,    83,    84,    85,    35,    92,   267,   130,   268,
     133,   148,   274,    96,    97,    98,    99,    81,    77,    78,
      79,   261,  -111,  -111,  -111,   136,   137,   162,   149,   218,
      87,    88,   138,   139,   140,   141,   142,   143,   144,   145,
     146,   171,   150,   184,   117,   106,   118,   119,   120,   121,
     122,   123,   113,   114,   185,   126,   168,   151,   127,   128,
     178,   257,   258,   259,    48,    49,   226,   227,   186,   187,
     188,   189,   190,   191,   192,   193,   194,   195,   196,   197,
     198,   199,     2,  -112,  -112,  -112,   232,   233,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    51,    52,    53,    54,    55,    56,    57,
      58,    59,    60,    61,    62,    63,    64,    65,    66,    67,
      68,    69,    70,    71,   202,   203,   204,   205,   206,   207,
     263,   264,   154,   156,   158,   163,   164,   159,   165,   180,
     166,   172,   173,   174,   175,   176,   209,   210,   211,   212,
     213,   214,   215,   216,   217,   220,   222,   223,   224,   234,
     235,   239,   219,   236,   237,   241,   242,   243,   244,   260,
     245,   246,   247,   248,   250,   249,   251,   253,   254,   262,
     266,   272,   275,   273,   277,   279,   278,   238,   152,   221,
     271,   153
};

static const yytype_uint16 yycheck[] =
{
       3,     7,    29,   181,    11,     7,     6,     6,     5,     6,
     229,     4,    42,    43,    50,    45,    11,     6,    33,    34,
      35,    36,    37,     7,    30,    45,    46,    47,    30,    48,
      49,     5,     6,   102,     6,    54,   255,    52,    53,    31,
      32,    78,    79,    80,    81,    37,    96,     3,    75,     5,
       5,    81,   271,    38,    39,    40,    41,    53,    65,    66,
      67,   239,   106,   107,   108,    62,    63,    97,     6,    69,
      65,    66,    69,    70,    71,    72,    73,    74,    75,    76,
      77,   108,     5,    57,    56,    42,    58,    59,    60,    61,
      62,    63,    43,    44,    68,    97,   102,     4,   100,   101,
      99,   106,   107,   108,    64,    65,   109,   110,    82,    83,
      84,    85,    86,    87,    88,    89,    90,    91,    92,    93,
      94,    95,     0,   106,   107,   108,     3,     4,     6,     7,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    49,    50,    51,    52,    53,    54,
     104,   105,     4,     4,     4,    35,    96,     6,    96,     4,
      96,    96,    96,    96,    96,    96,    96,    96,    96,    96,
      96,    96,    96,    96,    96,     5,     4,     4,     3,     3,
       3,    98,   103,     4,     4,     4,     3,     3,     3,    98,
       4,     3,     3,     3,    96,     4,     4,     3,   111,     3,
      98,   105,     3,   104,   104,     3,   105,   178,    86,   161,
     265,    87
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,   113,     0,     6,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,   114,   115,
     160,     4,   119,    31,    32,    37,   128,    33,    34,    35,
      36,    37,    52,    53,   133,    50,   138,   135,    64,    65,
     140,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,   150,     6,   151,   116,    11,    65,    66,    67,
     156,    53,    78,    79,    80,    81,    11,    65,    66,   160,
     102,   124,    96,   129,   129,   129,    38,    39,    40,    41,
     134,   129,    45,    46,    47,   139,    42,   136,   142,    48,
      49,    54,   137,    43,    44,   141,     6,    56,    58,    59,
      60,    61,    62,    63,   152,    30,    97,   100,   101,   117,
     160,   153,   158,     5,     5,     6,    62,    63,    69,    70,
      71,    72,    73,    74,    75,    76,    77,   157,   129,     6,
       5,     4,   153,   158,     4,   126,     4,   132,     4,     6,
     130,   131,   129,    35,    96,    96,    96,    30,   102,   143,
     144,   160,    96,    96,    96,    96,    96,     6,    99,   118,
       4,   120,     5,     6,    57,    68,    82,    83,    84,    85,
      86,    87,    88,    89,    90,    91,    92,    93,    94,    95,
     154,   155,    49,    50,    51,    52,    53,    54,   159,    96,
      96,    96,    96,    96,    96,    96,    96,    96,    69,   103,
       5,   132,     4,     4,     3,     3,   109,   110,   145,   146,
     147,   148,     3,     4,     3,     3,     4,     4,   119,    98,
     124,     4,     3,     3,     3,     4,     3,     3,     3,     4,
      96,     4,   127,     3,   111,   149,   145,   106,   107,   108,
      98,   124,     3,   104,   105,   145,    98,     3,     5,   121,
     125,   149,   105,   104,   145,     3,   122,   104,   105,     3,
     123
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


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
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

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
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}

/* Prevent warnings from -Wmissing-prototypes.  */
#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */


/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*-------------------------.
| yyparse or yypush_parse.  |
`-------------------------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{


    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks thru separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yytoken = 0;
  yyss = yyssa;
  yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */
  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

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

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 6:

/* Line 1455 of yacc.c  */
#line 114 "polypars"
    {  putchar('\n');
              putchar('\n');
              macrof = strdup(yylval.string);
              addmacro(polyres, con_b, con_f, no, macrof, phiv, psiv, omev);
              wline('-', MAXLEN);
           ;}
    break;

  case 7:

/* Line 1455 of yacc.c  */
#line 121 "polypars"
    {  putchar('\n'); 
              putchar('\n'); 
              tertiary(resno);
              Hydrolysis(da);
              BuildInfo();
              wline('-', MAXLEN);
           ;}
    break;

  case 8:

/* Line 1455 of yacc.c  */
#line 129 "polypars"
    {  putchar('\n');
              putchar('\n');
              calc_dist(Aijkl[0], Aijkl[1]);
              inx = 0;
              wline('-', MAXLEN);
           ;}
    break;

  case 9:

/* Line 1455 of yacc.c  */
#line 136 "polypars"
    {  putchar('\n');
              putchar('\n');
              calc_angle(Aijkl[0], Aijkl[1], Aijkl[2]);
              inx = 0;
              wline('-', MAXLEN);
           ;}
    break;

  case 10:

/* Line 1455 of yacc.c  */
#line 143 "polypars"
    {  putchar('\n');
              putchar('\n');
              calc_tors(Aijkl[0], Aijkl[1], Aijkl[2], Aijkl[3]);
              inx = 0;
              wline('-', MAXLEN);
           ;}
    break;

  case 11:

/* Line 1455 of yacc.c  */
#line 150 "polypars"
    {  putchar('\n');
              putchar('\n');
              dist_tabel();
              wline('-', MAXLEN);
           ;}
    break;

  case 12:

/* Line 1455 of yacc.c  */
#line 156 "polypars"
    {  putchar('\n');
              putchar('\n');
              calc_mass();
              wline('-', MAXLEN);
           ;}
    break;

  case 13:

/* Line 1455 of yacc.c  */
#line 162 "polypars"
    {  putchar('\n');
              putchar('\n');
              calc_rgyr();
              wline('-', MAXLEN);
           ;}
    break;

  case 14:

/* Line 1455 of yacc.c  */
#line 168 "polypars"
    {  putchar('\n');
              putchar('\n');
              calc_vol();
              wline('-', MAXLEN);
           ;}
    break;

  case 15:

/* Line 1455 of yacc.c  */
#line 174 "polypars"
    {  inx = 0;
           ;}
    break;

  case 16:

/* Line 1455 of yacc.c  */
#line 177 "polypars"
    {  putchar('\n');
              putchar('\n');
              add_bond(Aijkl[0],Aijkl[1]);
              calc_dist(Aijkl[0], Aijkl[1]);
              inx = 0;
              wline('-', MAXLEN);
           ;}
    break;

  case 17:

/* Line 1455 of yacc.c  */
#line 185 "polypars"
    {  putchar('\n');
              putchar('\n');
              pdate(tid);
              wline('-', MAXLEN);
           ;}
    break;

  case 18:

/* Line 1455 of yacc.c  */
#line 191 "polypars"
    {  putchar('\n');
              putchar('\n');
              doubleh(nstr,nreps,shift);
              wline('-', MAXLEN);  
           ;}
    break;

  case 19:

/* Line 1455 of yacc.c  */
#line 197 "polypars"
    {  printf("\n\n\tbye - bye\n\n");
              return (0); 
           ;}
    break;

  case 20:

/* Line 1455 of yacc.c  */
#line 201 "polypars"
    {  putchar('\n');
              putchar('\n');
              printener(deriv); 
              deriv = 0;
              wline('-', MAXLEN);
           ;}
    break;

  case 21:

/* Line 1455 of yacc.c  */
#line 208 "polypars"
    {  putchar('\n');
              putchar('\n');
              Generate(pargenopt); 
              wline('-', MAXLEN);
           ;}
    break;

  case 22:

/* Line 1455 of yacc.c  */
#line 214 "polypars"
    {  putchar('\n');
              wline('-', MAXLEN);  
           ;}
    break;

  case 23:

/* Line 1455 of yacc.c  */
#line 218 "polypars"
    {  putchar('\n');
              putchar('\n');
              calc_helix(optimize, nfold);
              wline('-', MAXLEN);  
           ;}
    break;

  case 24:

/* Line 1455 of yacc.c  */
#line 224 "polypars"
    {  putchar('\n');
              putchar('\n');
              printf("Initializing data structures......\n");
              Initialize();
              wline('-', MAXLEN);
           ;}
    break;

  case 25:

/* Line 1455 of yacc.c  */
#line 231 "polypars"
    {  putchar('\n');
              putchar('\n');
              phelp(yylval.cmd);
              wline('-', MAXLEN);
           ;}
    break;

  case 26:

/* Line 1455 of yacc.c  */
#line 237 "polypars"
    {  putchar('\n');
              putchar('\n');
              printf("Label: '%s' in line %d\n", yylval.string, lineno);
              wline('-', MAXLEN);
           ;}
    break;

  case 27:

/* Line 1455 of yacc.c  */
#line 243 "polypars"
    {  putchar('\n');
              putchar('\n');
              Minimize(parspace,maxfun,maxit,enetol,grdtol,silent);
              parspace = 0;
              wline('-', MAXLEN);
           ;}
    break;

  case 28:

/* Line 1455 of yacc.c  */
#line 250 "polypars"
    {  putchar('\n');
              putchar('\n');
              phiv = yylval.real;
              Rotate(Aijkl[0], Aijkl[1], Aijkl[2], Aijkl[3], phiv);
              inx = 0;
              wline('-', MAXLEN);
           ;}
    break;

  case 29:

/* Line 1455 of yacc.c  */
#line 258 "polypars"
    {  putchar('\n'); 
              putchar('\n'); 
              res_list(resno);
              wline('-', MAXLEN);
           ;}
    break;

  case 30:

/* Line 1455 of yacc.c  */
#line 264 "polypars"
    {  putchar('\n');
              putchar('\n');
              fname = strdup(yylval.string);
              ReadCharges(fname);
              free(fname);
              wline('-', MAXLEN);
           ;}
    break;

  case 31:

/* Line 1455 of yacc.c  */
#line 272 "polypars"
    {  putchar('\n');
              putchar('\n');
              ReadCoord(parformat,fname,reference);
              free(fname);
              wline('-', MAXLEN);
           ;}
    break;

  case 32:

/* Line 1455 of yacc.c  */
#line 279 "polypars"
    {  putchar('\n');
              putchar('\n');
              printf("Reading internal coordinates......\n");
              wline('-', MAXLEN);
           ;}
    break;

  case 33:

/* Line 1455 of yacc.c  */
#line 285 "polypars"
    {  putchar('\n');
              putchar('\n');
              ReadParam();
              wline('-', MAXLEN);
           ;}
    break;

  case 34:

/* Line 1455 of yacc.c  */
#line 291 "polypars"
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
           ;}
    break;

  case 35:

/* Line 1455 of yacc.c  */
#line 304 "polypars"
    {  putchar('\n');
              putchar('\n');
              debug = TRUE;
              printf("Debugging mode on\n");
              wline('-', MAXLEN);
           ;}
    break;

  case 36:

/* Line 1455 of yacc.c  */
#line 311 "polypars"
    {  putchar('\n');
              putchar('\n');
              M.id = strdup(yylval.string);
              printf("ID: >%s<\n", M.id);
              wline('-', MAXLEN);
           ;}
    break;

  case 37:

/* Line 1455 of yacc.c  */
#line 318 "polypars"
    {  putchar('\n');
              putchar('\n');
              seed = -(yylval.integer);
              printf("The random number generator has been initialized with");
              printf(" a new seed (%d)\n", yylval.integer);
              wline('-', MAXLEN);
           ;}
    break;

  case 38:

/* Line 1455 of yacc.c  */
#line 326 "polypars"
    {  putchar('\n');
              putchar('\n');
              M.titel = strdup(yylval.string);
              printf("TITLE: %s\n", M.titel);
              wline('-', MAXLEN);
           ;}
    break;

  case 39:

/* Line 1455 of yacc.c  */
#line 333 "polypars"
    {
           ;}
    break;

  case 40:

/* Line 1455 of yacc.c  */
#line 336 "polypars"
    {  putchar('\n');
              putchar('\n');
              ptime(tid);
              wline('-', MAXLEN);
           ;}
    break;

  case 41:

/* Line 1455 of yacc.c  */
#line 342 "polypars"
    {  putchar('\n');
              putchar('\n');
              WritCoord(parformat, fname, descarte); 
              descarte = TRUE;
              free(fname);
              wline('-', MAXLEN);
           ;}
    break;

  case 42:

/* Line 1455 of yacc.c  */
#line 350 "polypars"
    {  putchar('\n');
              putchar('\n');
              list_intcoor(intopt); 
              wline('-', MAXLEN);
           ;}
    break;

  case 43:

/* Line 1455 of yacc.c  */
#line 356 "polypars"
    {  putchar('\n');
              putchar('\n');
              printf("\nWriting parameters ........\n"); 
              wline('-', MAXLEN);
           ;}
    break;

  case 44:

/* Line 1455 of yacc.c  */
#line 362 "polypars"
    {  putchar('\n');
              putchar('\n');
              yyerror("unknown command:"); 
              wline('-', MAXLEN);
           ;}
    break;

  case 48:

/* Line 1455 of yacc.c  */
#line 373 "polypars"
    { resno++;
             expand(resno,fname,con_b,con_f,thev,phiv,psiv,omev,resno+1,status);
             free(fname);
           ;}
    break;

  case 49:

/* Line 1455 of yacc.c  */
#line 378 "polypars"
    {  resno++;
              expand(resno,fname,0,0,999.9,999.9,999.9,999.9,0,status);
              free(fname);
           ;}
    break;

  case 50:

/* Line 1455 of yacc.c  */
#line 383 "polypars"
    {  oldstatus=status;
              status++;
              /* if we are in a repeat sequence then save fconnect */
              if ((repno > 0) AND (reptab[repno].fconnect==0))
                 reptab[repno].fconnect = resno;
           ;}
    break;

  case 51:

/* Line 1455 of yacc.c  */
#line 390 "polypars"
    {   
              modify_con(resno,polyres,status);
              status=oldstatus;
           ;}
    break;

  case 52:

/* Line 1455 of yacc.c  */
#line 395 "polypars"
    {  repno++;
              reptab[repno].fconnect=0;
              reptab[repno].last=0;
              reptab[repno].rep=0;
              reptab[repno].first = resno+1;
           ;}
    break;

  case 53:

/* Line 1455 of yacc.c  */
#line 402 "polypars"
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
           ;}
    break;

  case 54:

/* Line 1455 of yacc.c  */
#line 428 "polypars"
    {  fname = strdup(yylval.string); ;}
    break;

  case 55:

/* Line 1455 of yacc.c  */
#line 431 "polypars"
    {  polyres = yylval.integer; 
           ;}
    break;

  case 56:

/* Line 1455 of yacc.c  */
#line 435 "polypars"
    {  reptab[repno].rep = yylval.integer; ;}
    break;

  case 57:

/* Line 1455 of yacc.c  */
#line 438 "polypars"
    {  phiv = 999.9;
              psiv = 999.9;
              omev = 999.9;
              phiv = yylval.real; ;}
    break;

  case 58:

/* Line 1455 of yacc.c  */
#line 444 "polypars"
    {  psiv = yylval.real; ;}
    break;

  case 59:

/* Line 1455 of yacc.c  */
#line 447 "polypars"
    {  omev = yylval.real; ;}
    break;

  case 64:

/* Line 1455 of yacc.c  */
#line 455 "polypars"
    {  phiv = 999.9;
              psiv = 999.9;
              omev = 999.9;
              mapfile = strdup(yylval.string); 
              getconval(mapfile,&phiv,&psiv,&omev);
              free(mapfile);
           ;}
    break;

  case 65:

/* Line 1455 of yacc.c  */
#line 464 "polypars"
    {  con_b = yylval.integer; ;}
    break;

  case 66:

/* Line 1455 of yacc.c  */
#line 467 "polypars"
    {  con_f = yylval.integer; ;}
    break;

  case 67:

/* Line 1455 of yacc.c  */
#line 470 "polypars"
    {  nsets  =  yylval.integer; ;}
    break;

  case 68:

/* Line 1455 of yacc.c  */
#line 472 "polypars"
    { ;}
    break;

  case 69:

/* Line 1455 of yacc.c  */
#line 474 "polypars"
    { ;}
    break;

  case 73:

/* Line 1455 of yacc.c  */
#line 481 "polypars"
    {  Aijkl[inx++] = getAno(no, lab); ;}
    break;

  case 74:

/* Line 1455 of yacc.c  */
#line 484 "polypars"
    {  strcpy(lab, yyval.string); ;}
    break;

  case 75:

/* Line 1455 of yacc.c  */
#line 487 "polypars"
    {  no = yyval.integer; ;}
    break;

  case 78:

/* Line 1455 of yacc.c  */
#line 493 "polypars"
    {  Coord_copy(); ;}
    break;

  case 79:

/* Line 1455 of yacc.c  */
#line 495 "polypars"
    {  putchar('\n');
              putchar('\n');
              Coord_orient(Aijkl[0],Aijkl[1],Aijkl[2]);
              wline('-', MAXLEN);
           ;}
    break;

  case 80:

/* Line 1455 of yacc.c  */
#line 501 "polypars"
    {  putchar('\n');
              putchar('\n');
              Coord_rms(FALSE);
              wline('-', MAXLEN);
           ;}
    break;

  case 81:

/* Line 1455 of yacc.c  */
#line 507 "polypars"
    {  putchar('\n');
              putchar('\n');
              Coord_rms(TRUE);
              wline('-', MAXLEN);
           ;}
    break;

  case 82:

/* Line 1455 of yacc.c  */
#line 513 "polypars"
    {  putchar('\n');
              putchar('\n');
              Coord_stat(); 
              wline('-', MAXLEN);
           ;}
    break;

  case 85:

/* Line 1455 of yacc.c  */
#line 523 "polypars"
    {  deriv = 1; ;}
    break;

  case 86:

/* Line 1455 of yacc.c  */
#line 525 "polypars"
    {  pargenopt = 1; ;}
    break;

  case 87:

/* Line 1455 of yacc.c  */
#line 527 "polypars"
    {  pargenopt = 0; ;}
    break;

  case 88:

/* Line 1455 of yacc.c  */
#line 529 "polypars"
    {  pargenopt = 1; ;}
    break;

  case 89:

/* Line 1455 of yacc.c  */
#line 531 "polypars"
    {  pargenopt = 11; ;}
    break;

  case 92:

/* Line 1455 of yacc.c  */
#line 537 "polypars"
    {  nstr = yylval.integer; ;}
    break;

  case 93:

/* Line 1455 of yacc.c  */
#line 539 "polypars"
    {  nreps = yylval.integer; ;}
    break;

  case 94:

/* Line 1455 of yacc.c  */
#line 541 "polypars"
    {  shift = yylval.real; ;}
    break;

  case 97:

/* Line 1455 of yacc.c  */
#line 547 "polypars"
    {  nfold = yylval.real; ;}
    break;

  case 98:

/* Line 1455 of yacc.c  */
#line 549 "polypars"
    {  nfold = (double) yylval.integer; ;}
    break;

  case 99:

/* Line 1455 of yacc.c  */
#line 551 "polypars"
    {  optimize = TRUE; ;}
    break;

  case 102:

/* Line 1455 of yacc.c  */
#line 557 "polypars"
    {  if (rnat == 0)
                 rnat = M.nat;
              symneighbor(rnat, TM); 
           ;}
    break;

  case 104:

/* Line 1455 of yacc.c  */
#line 564 "polypars"
    {  initmat(TM, 0.0); 
              TM[3][3] = 1.0;
              symno = 0;
           ;}
    break;

  case 107:

/* Line 1455 of yacc.c  */
#line 573 "polypars"
    {  TM[symno][0] = xmult; ;}
    break;

  case 108:

/* Line 1455 of yacc.c  */
#line 575 "polypars"
    {  TM[symno][1] = xmult; ;}
    break;

  case 109:

/* Line 1455 of yacc.c  */
#line 577 "polypars"
    {  TM[symno][2] = xmult; ;}
    break;

  case 110:

/* Line 1455 of yacc.c  */
#line 579 "polypars"
    {  TM[symno][3] = xtransl; ;}
    break;

  case 111:

/* Line 1455 of yacc.c  */
#line 582 "polypars"
    {  xmult = 1.0; ;}
    break;

  case 112:

/* Line 1455 of yacc.c  */
#line 584 "polypars"
    {  xmult = yylval.real; ;}
    break;

  case 113:

/* Line 1455 of yacc.c  */
#line 586 "polypars"
    {  xmult = 1.0; ;}
    break;

  case 114:

/* Line 1455 of yacc.c  */
#line 588 "polypars"
    {  xmult = -1.0; ;}
    break;

  case 115:

/* Line 1455 of yacc.c  */
#line 591 "polypars"
    {  xtransl = yylval.real; ;}
    break;

  case 116:

/* Line 1455 of yacc.c  */
#line 593 "polypars"
    {  xtransl = yylval.real; ;}
    break;

  case 117:

/* Line 1455 of yacc.c  */
#line 596 "polypars"
    {  symno++; ;}
    break;

  case 119:

/* Line 1455 of yacc.c  */
#line 600 "polypars"
    { ;}
    break;

  case 120:

/* Line 1455 of yacc.c  */
#line 602 "polypars"
    { ;}
    break;

  case 121:

/* Line 1455 of yacc.c  */
#line 604 "polypars"
    { ;}
    break;

  case 122:

/* Line 1455 of yacc.c  */
#line 606 "polypars"
    { ;}
    break;

  case 123:

/* Line 1455 of yacc.c  */
#line 608 "polypars"
    { ;}
    break;

  case 124:

/* Line 1455 of yacc.c  */
#line 610 "polypars"
    { ;}
    break;

  case 125:

/* Line 1455 of yacc.c  */
#line 612 "polypars"
    { ;}
    break;

  case 126:

/* Line 1455 of yacc.c  */
#line 614 "polypars"
    { ;}
    break;

  case 127:

/* Line 1455 of yacc.c  */
#line 616 "polypars"
    { ;}
    break;

  case 128:

/* Line 1455 of yacc.c  */
#line 618 "polypars"
    { ;}
    break;

  case 129:

/* Line 1455 of yacc.c  */
#line 620 "polypars"
    { ;}
    break;

  case 130:

/* Line 1455 of yacc.c  */
#line 622 "polypars"
    { ;}
    break;

  case 131:

/* Line 1455 of yacc.c  */
#line 624 "polypars"
    { ;}
    break;

  case 132:

/* Line 1455 of yacc.c  */
#line 626 "polypars"
    { ;}
    break;

  case 133:

/* Line 1455 of yacc.c  */
#line 628 "polypars"
    { ;}
    break;

  case 134:

/* Line 1455 of yacc.c  */
#line 630 "polypars"
    { ;}
    break;

  case 135:

/* Line 1455 of yacc.c  */
#line 632 "polypars"
    { ;}
    break;

  case 136:

/* Line 1455 of yacc.c  */
#line 634 "polypars"
    { ;}
    break;

  case 137:

/* Line 1455 of yacc.c  */
#line 636 "polypars"
    { ;}
    break;

  case 138:

/* Line 1455 of yacc.c  */
#line 638 "polypars"
    { ;}
    break;

  case 139:

/* Line 1455 of yacc.c  */
#line 640 "polypars"
    { ;}
    break;

  case 142:

/* Line 1455 of yacc.c  */
#line 646 "polypars"
    {  parspace = 1; ;}
    break;

  case 143:

/* Line 1455 of yacc.c  */
#line 648 "polypars"
    {  grdtol =  yylval.real; ;}
    break;

  case 144:

/* Line 1455 of yacc.c  */
#line 650 "polypars"
    {  enetol =  yylval.real; ;}
    break;

  case 145:

/* Line 1455 of yacc.c  */
#line 652 "polypars"
    {  maxit  =  yylval.integer; ;}
    break;

  case 146:

/* Line 1455 of yacc.c  */
#line 654 "polypars"
    {  maxfun =  yylval.integer; ;}
    break;

  case 147:

/* Line 1455 of yacc.c  */
#line 656 "polypars"
    {  silent = FALSE; ;}
    break;

  case 148:

/* Line 1455 of yacc.c  */
#line 658 "polypars"
    {  silent = TRUE; ;}
    break;

  case 149:

/* Line 1455 of yacc.c  */
#line 660 "polypars"
    { yyerror("unknown subcommand: "); ;}
    break;

  case 153:

/* Line 1455 of yacc.c  */
#line 667 "polypars"
    {  descarte = FALSE; ;}
    break;

  case 154:

/* Line 1455 of yacc.c  */
#line 669 "polypars"
    {  reference = TRUE; ;}
    break;

  case 155:

/* Line 1455 of yacc.c  */
#line 671 "polypars"
    {  fname = strdup(yylval.string); ;}
    break;

  case 156:

/* Line 1455 of yacc.c  */
#line 673 "polypars"
    { yyerror("unknown subcommand: "); ;}
    break;

  case 157:

/* Line 1455 of yacc.c  */
#line 676 "polypars"
    { parformat = 0; ;}
    break;

  case 158:

/* Line 1455 of yacc.c  */
#line 678 "polypars"
    { parformat = 1; ;}
    break;

  case 159:

/* Line 1455 of yacc.c  */
#line 680 "polypars"
    { parformat = 2; ;}
    break;

  case 160:

/* Line 1455 of yacc.c  */
#line 682 "polypars"
    { parformat = 3; ;}
    break;

  case 161:

/* Line 1455 of yacc.c  */
#line 684 "polypars"
    { parformat = 4; ;}
    break;

  case 162:

/* Line 1455 of yacc.c  */
#line 686 "polypars"
    { parformat = 5; ;}
    break;

  case 163:

/* Line 1455 of yacc.c  */
#line 688 "polypars"
    { parformat = 6; ;}
    break;

  case 164:

/* Line 1455 of yacc.c  */
#line 690 "polypars"
    { parformat = 7; ;}
    break;

  case 165:

/* Line 1455 of yacc.c  */
#line 692 "polypars"
    { parformat = 8; ;}
    break;

  case 166:

/* Line 1455 of yacc.c  */
#line 694 "polypars"
    { parformat = 9; ;}
    break;

  case 167:

/* Line 1455 of yacc.c  */
#line 696 "polypars"
    { parformat = 10; ;}
    break;

  case 168:

/* Line 1455 of yacc.c  */
#line 698 "polypars"
    { parformat = 11; ;}
    break;

  case 169:

/* Line 1455 of yacc.c  */
#line 700 "polypars"
    { parformat = 12; ;}
    break;

  case 170:

/* Line 1455 of yacc.c  */
#line 702 "polypars"
    { parformat = 13; ;}
    break;

  case 173:

/* Line 1455 of yacc.c  */
#line 708 "polypars"
    {  phino = yylval.integer; ;}
    break;

  case 174:

/* Line 1455 of yacc.c  */
#line 710 "polypars"
    {  phib  = yylval.real; ;}
    break;

  case 175:

/* Line 1455 of yacc.c  */
#line 712 "polypars"
    {  phie  = yylval.real; ;}
    break;

  case 176:

/* Line 1455 of yacc.c  */
#line 714 "polypars"
    {  phiinc= yylval.real; ;}
    break;

  case 177:

/* Line 1455 of yacc.c  */
#line 716 "polypars"
    {  psino = yylval.integer; ;}
    break;

  case 178:

/* Line 1455 of yacc.c  */
#line 718 "polypars"
    {  psib  = yylval.real; ;}
    break;

  case 179:

/* Line 1455 of yacc.c  */
#line 720 "polypars"
    {  psie  = yylval.real; ;}
    break;

  case 180:

/* Line 1455 of yacc.c  */
#line 722 "polypars"
    {  psiinc= yylval.real; ;}
    break;

  case 181:

/* Line 1455 of yacc.c  */
#line 724 "polypars"
    {  recopt = yylval.integer; ;}
    break;

  case 182:

/* Line 1455 of yacc.c  */
#line 726 "polypars"
    {  silent = FALSE; ;}
    break;

  case 183:

/* Line 1455 of yacc.c  */
#line 728 "polypars"
    {  silent = TRUE; ;}
    break;

  case 184:

/* Line 1455 of yacc.c  */
#line 730 "polypars"
    {  fname = strdup(yylval.string); ;}
    break;

  case 185:

/* Line 1455 of yacc.c  */
#line 732 "polypars"
    { yyerror("unknown subcommand: "); ;}
    break;

  case 186:

/* Line 1455 of yacc.c  */
#line 734 "polypars"
    { intopt = 1; ;}
    break;

  case 188:

/* Line 1455 of yacc.c  */
#line 738 "polypars"
    { intopt = 2; ;}
    break;

  case 189:

/* Line 1455 of yacc.c  */
#line 740 "polypars"
    { intopt = 1; ;}
    break;

  case 190:

/* Line 1455 of yacc.c  */
#line 742 "polypars"
    { intopt = 1; ;}
    break;

  case 191:

/* Line 1455 of yacc.c  */
#line 744 "polypars"
    { intopt = 1; ;}
    break;

  case 192:

/* Line 1455 of yacc.c  */
#line 746 "polypars"
    { intopt = 1; ;}
    break;

  case 193:

/* Line 1455 of yacc.c  */
#line 748 "polypars"
    { intopt = 0; ;}
    break;

  case 194:

/* Line 1455 of yacc.c  */
#line 751 "polypars"
    {  lineno++;
              if (lineno > 10000) lineno=1;
              printf("\nPOLYS:%d > ", lineno);
           ;}
    break;



/* Line 1455 of yacc.c  */
#line 3136 "polypars.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

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

#if !defined(yyoverflow) || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}



/* Line 1675 of yacc.c  */
#line 756 "polypars"



