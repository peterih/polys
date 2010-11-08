
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton interface for Bison's Yacc-like parsers in C
   
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

/* Line 1676 of yacc.c  */
#line 82 "polypars"
 double   real;
         int      integer;
         char     *string;
         int      cmd;
       


/* Line 1676 of yacc.c  */
#line 155 "polypars.tab.h"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE yylval;


