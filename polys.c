/******************************************************************POLYS
** polys.c                                                            **
**                                                                    **
** This module contains the subroutine "main" to programme polys:     **
**                                                                    **
** THE POLYSaccharide BUILDER                                         **
**                                                                    **
** S.B. Engelsen, S. Cros, W. Mackie & S. Pérez,                      **
**  A molecular builder for carbohydrates.                            **
**  Application to polysaccharides and complex carbohydrates,         ** 
**  Biopolymers (1996), 39, 417-433.                                  **
**                                                                    **
** Written by Soren Balling Engelsen, INRA-93/94                      **
** Revised by Soren Balling Engelsen, UCPH-97/12                      **
***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "polys.h"
#include "vars.h"

/* Print welcome header on output device */
void
header(void)
{
   printf("\n\n");
   printf("\tPPPP   OOO  L    Y   Y  SSS\n");
   printf("\tP   P O   O L     Y Y  S\n");
   printf("\tPPPP  O   O L      Y    SSS\n");
   printf("\tP     O   O L      Y       S\n");
   printf("\tP      OOO  LLLL   Y   SSSS\n");
   printf("\n");
   printf("\tTHE POLYSaccharide BUILDER\n\n");
   printf("\tSoren Balling Engelsen and Serge Perez\n\n");
   printf("\tNantes - Copenhagen - Grenoble\n\n");
   printf("\tPOLYS VERSION 2.0 (2012)\n\n");
   printf("\tCode Contributors:\n\n");
   printf("\tSoren Balling Engelsen * ");
   printf("Milou Kouwijzer * ");
   printf("Peter Ibsen Hansen\n");
}

/* Print welcome header on output device */
main(int argc, char **argv)
{
   int        yy, mm, dd, hh, min, ss;
   double     tforce;
   
   /* initialize time */
   tid = 0.0;
   tid = etime(tid);
 
   ProgName = argv[0];
   UserName = getenv("USER");

   /* Initialize datastructures and write header */
   Initialize();
   wline('#', MAXLEN);
   header();
   wline('#', MAXLEN);
   printf("\nPOLYS:%d > ", lineno);

   /* invoke polys command line parser */
   yyparse();

   /* write session tail */
   dtime(&yy,&mm,&dd,&hh,&min,&ss);
   printf("\nProgram %s ending >%s< session for user: %s\n\n",
                   ProgName,   M.id,             UserName);
   printf("Session ended at %d:%02d on %d-%d-%d\n\n", hh, min, dd, mm, yy);
   ttot = etime(0.0);
   printf("Elapsed CPU time:     Function     Seconds         %%        %%\n");
   printf("----------------------------------------------------------------\n");
   printf("%30s%12.2lf%12.2lf\n", "Builder", tbuild, tbuild/ttot*100.0);
   printf("%30s%12.2lf%12.2lf\n", "Internal", tint, tint/ttot*100.0);
   printf("%30s%12.2lf%12.2lf\n", "Param", tpara, tpara/ttot*100.0);
   tforce = tnb + tbo + tv + tt;
   printf("%30s%12.2lf%12.2lf%9.2lf\n", "Bond", tbo, tbo/ttot*100.0, tbo/tforce*100.0);
   printf("%30s%12.2lf%12.2lf%9.2lf\n", "Valence", tv, tv/ttot*100.0, tv/tforce*100.0);
   printf("%30s%12.2lf%12.2lf%9.2lf\n", "Torsion", tt, tt/ttot*100.0, tt/tforce*100.0);
   printf("%30s%12.2lf%12.2lf%9.2lf\n", "Nonbond", tnb, tnb/ttot*100.0, tnb/tforce*100.0);
   printf("%30s%12.2lf%12.2lf\n", "Rotation", trot, trot/ttot*100.0);
   tother = ttot - (tbuild + tint + tpara + tnb + tbo + tv + tt + trot); 
   printf("%30s%12.2lf%12.2lf\n", "Other", tother, tother/ttot*100.0);
   if (ttest > 0.0)
   printf("%30s%12.2lf%12.2lf\n", "Test", ttest, ttest/ttot*100.0);
   printf("----------------------------------------------------------------\n");
   printf("%30s%12.2lf%12.2lf%9.2lf\n\n", "Total", ttot, 100.0, 100.0);
   wline('#', MAXLEN);
   printf("\n\n");

   return (0);
}

/* Print syntax error from the parser */
void
yyerror(char *s)
{
   fprintf(stderr, "\n%s: %s\n", ProgName, s);
   fprintf(stderr, " line %d\n\n", lineno);
}
