#ifndef BISON_POLYPARS_TAB_H
# define BISON_POLYPARS_TAB_H

#ifndef YYSTYPE
typedef union { double   real;
         int      integer;
         char     *string;
         int      cmd;
       } yystype;
# define YYSTYPE yystype
# define YYSTYPE_IS_TRIVIAL 1
#endif
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


extern YYSTYPE yylval;

#endif /* not BISON_POLYPARS_TAB_H */
