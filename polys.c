/**********************************************************************
POLYS MAIN MODULE <polys.c>

   This module contains the subroutine "main" to program polys.

   Written by Soren Balling Engelsen, INRA-93.
   Last revision: S.B.E. 1997
**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "polys.h"
#include "vars.h"


/**********************************************************************
	Print welcome header on output device.
**********************************************************************/
void header(void)

{  printf("\n\n\n\n");
   printf("\tPPPPPPP      OOOOO     LL      YY      YY   SSSSSS\n");
   printf("\tPP   PPP   OOO   OOO   LL       YY    YY   SSS\n");
   printf("\tPP   PPP  OO       OO  LL        YY  YY    SS\n");
   printf("\tPPPPPPP   OO       OO  LL         YYYY      SSSSSS\n");
   printf("\tPP        OO       OO  LL          YY            SS\n");
   printf("\tPP         OOO   OOO   LL          YY           SSS\n");
   printf("\tPP           OOOOO     LLLLLLL     YY       SSSSSS\n");
   printf("\n\n\n\n\n");
   printf("\tTHE POLYSaccharide BUILDER (1993-95) - INRA Nantes\n\n");
   printf("\t\tSoren Balling Engelsen and Serge Perez\n\n\n");
   printf("\t\tKVL                        CERMAV\n\n\n");
   printf("\t\tPOLYS VERSION 1.4 (2004)\n\n");
   printf("\t\tCode Contributors\n\n");
   printf("\t\tSoren Balling Engelsen\n\n");
   printf("\t\tMilou Kouwijzer\n\n");
   printf("\t\tPeter Ibsen Hansen\n\n");

} /* End of header */


/**********************************************************************
	Print welcome header on output device.
**********************************************************************/
main(int argc, char **argv)
{  int        yy, mm, dd, hh, min, ss;
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

}  /* End of main */


/**********************************************************************
   Print syntax error from the parser.
**********************************************************************/
void yyerror(char *s)

{  fprintf(stderr, "\n%s: %s\n", ProgName, s);
   fprintf(stderr, " line %d\n\n", lineno);
   exit(0);

} /* End of yyerror */


/* End of file */
