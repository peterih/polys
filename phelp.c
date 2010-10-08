/**********************************************************************
  "phelp.c"

  This module contains the source and text for help functions.

  Written by Soren Balling Engelsen, INRA-95.
**********************************************************************/
#include <string.h>
#include "polypars.tab.h"
 

/**********************************************************************
<cmd> ATTACH BUILD CALC COORD CREATE DATE DOUBLEHELIX END ENER GENERATE INIT
<cmd> HELIX HELP LABEL MINI PRIMARY READ RECORD MODIFY SET TIME WRITE
**********************************************************************/
void phelp(int ycom)

{  register int i;

   printf("POLYS HELP INFORMATION\n");

   switch (ycom)
   {  case    -1 :
      case ATTACH : 
       printf("\nThe ATTACH command\n");
       printf("SYNTAX:\tATTACH polyres=<integer> [connexion] macrores=<integer> <macro filename> \n");
       printf("FUNCTION:\tAttachment of a macro structure from a file\n");
       printf("FUNCTION:\tPrior to the execution of this command a primary\n");
       printf("FUNCTION:\tsequence must have been build using the PRIMARY\n");
       printf("FUNCTION:\tcommand\n");
       if (ycom != -1) break;
      case BUILD : 
       printf("\nThe BUILD command\n");
       printf("SYNTAX:\tBUILD\n");
       printf("FUNCTION:\tInvokes the polysaccharide building option.\n");
       printf("FUNCTION:\tPrior to the execution of this command a primary\n");
       printf("FUNCTION:\tsequence must have been build using the PRIMARY\n");
       printf("FUNCTION:\tcommand\n");
       if (ycom != -1) break;
      case CALC  :
       printf("\nThe CALC command\n");
       printf("SYNTAX:\tCALC [dist|angle|torsion] [atom]+\n");
       printf("SYNTAX:\tCALC [dtab|mass|rgyr|volume]\n");
       printf("FUNCTION:\tCalculates structural parameters:\n\n");
       printf("SUBCMD:\tdist [atom]+\n");
       printf("FUNCTION:\tInteratomic distance\n");
       printf("EX:\tCALC dist O_1 1 O_4 2\n\n");
       printf("SUBCMD: angle [atom]+\n");
       printf("FUNCTION:\tAngle defined by three atoms\n");
       printf("EX:\tCALC angle C_1 1 O_1 1 C_4 2\n\n");
       printf("SUBCMD: torsion [atom]+\n");
       printf("FUNCTION:\tTorsion defined by four atoms\n");
       printf("EX:\tCALC torsion O_5 1 C_1 1 O_1 1 C_4 2\n\n");
       printf("SUBCMD: dtab\n");
       printf("FUNCTION:\tInteratomic distance tabel\n");
       printf("EX:\tCALC dtab\n\n");
       printf("SUBCMD: mass\n");
       printf("FUNCTION:\tMolecular mass\n");
       printf("EX:\tCALC mass\n\n");
       printf("SUBCMD: rgyr\n");
       printf("FUNCTION:\tRadius of Gyration\n");
       printf("EX:\tCALC rgyr\n\n");
       printf("SUBCMD: volume\n");
       printf("FUNCTION:\tMolecular Volume\n");
       printf("EX:\tCALC volume\n");
       if (ycom != -1) break;
      case COORD :
       printf("\nThe COORD command\n");
       printf("SYNTAX:\tCOORD [copy|rms|rms mass|stat]\n");
       printf("SYNTAX:\tCOORD orient [atom]+\n");
       printf("FUNCTION:\tPerforms coordinate calculations and manipulations:\n\n");
       printf("SUBCMD:\tcopy\n");
       printf("FUNCTION:\tCopying coordinates to reference set\n");
       printf("EX:\tCOORD copy\n\n");
       printf("SUBCMD:\trms\n");
       printf("SUBCMD:\trms mass\n");
       printf("SUBCMD:\tstat\n");
       printf("SUBCMD:\torient [atom]+\n");
       if (ycom != -1) break;
      case CREATE:
       printf("\nThe CREATE command\n");
       printf("SYNTAX:\tCREATE [bond|angle|torsion] [atom]+\n");
       printf("FUNCTION:\tCreates new entries in the internal coordinate\n");
       printf("FUNCTION:\ttables:\n\n");
       printf("SUBCMD:\tbond [atom]+\n");
       printf("FUNCTION:\tNew bond definition\n");
       printf("EX:\tCREATE bond O_1 1 HO_1 1\n\n");
       printf("SUBCMD: angle [atom]+\n");
       printf("FUNCTION:\tNew valence angle definition\n");
       printf("EX:\tCREATE angle C_1 1 O_1 1 HO_1 1\n\n");
       printf("SUBCMD: torsion [atom]+\n");
       printf("FUNCTION:\tNew torsion definition\n");
       printf("EX:\tCREATE torsion O_5 1 C_1 1 O_1 1 HO_1 1\n\n");
       if (ycom != -1) break;
      case DATE  :
       printf("\nThe DATE command\n");
       printf("SYNTAX:\tDATE\n");
       printf("FUNCTION:\tWrites time and date to standard output:\n\n");
       if (ycom != -1) break;
      case DOUBLEHELIX:
       printf("\nThe DOUBLEHELIX command\n");
       printf("SYNTAX:\tDOUBLEHELIX nstr=<integer> [nreps=<integer>]\n");
       printf("SYNTAX:\t            [shift=<double>]\n");
       printf("FUNCTION:\t    \n\n");
       if (ycom != -1) break;
      case END:
       printf("\nThe END command\n");
       printf("SYNTAX:\tEND\n");
       printf("FUNCTION:\tEnds the POLYS session.\n\n");
       if (ycom != -1) break;
      case ENER:
       printf("\nThe ENER command\n");
       printf("SYNTAX:\tENER [deriv]\n");
       printf("FUNCTION:\tCalculates the potential energy of the molecular\n");
       printf("FUNCTION:\tsystem using the force field loaded with the \n");
       printf("FUNCTION:\tREAD param command\n\n");
       if (ycom != -1) break;
      case GENERATE:
       printf("\nThe GENERATE command\n");
       printf("SYNTAX:\tGENERATE [one|all|rotbond]\n");
       printf("FUNCTION:\t    \n\n");
       if (ycom != -1) break;
      case INIT:
       printf("\nThe INIT command\n");
       printf("SYNTAX:\tINIT\n");
       printf("FUNCTION:\t    \n\n");
       if (ycom != -1) break;
      case HELIX:
       printf("\nThe HELIX command\n");
       printf("SYNTAX:\tHELIX [optimize] [nfold=<integer>|<real>]\n");
       printf("FUNCTION:\t    \n\n");
       if (ycom != -1) break;
      case HELP:
       printf("\nThe HELP command\n");
       printf("SYNTAX:\tHELP [command]\n");
       printf("FUNCTION:\t    \n\n");
       if (ycom != -1) break;
      case LABEL:
       printf("\nThe LABEL command\n");
       printf("SYNTAX:\tLABEL <name>]\n");
       printf("FUNCTION:\t    \n\n");
       if (ycom != -1) break;
      case MINI:
       printf("\nThe MINI command\n");
       printf("SYNTAX:\tMINI [dihed] [noise|silent] [grdtol <value>]\n");
       printf("SYNTAX:\t     [enetol <value>] [iter <value>]\n");
       printf("FUNCTION:\t    \n\n");
       if (ycom != -1) break;
      case MODIFY:
       printf("\nThe MODIFY command\n");
       printf("SYNTAX:\tMODIFY [angle|torsion] [atom]+ [phi=<value>]\n");
       printf("FUNCTION:\t    \n\n");
       if (ycom != -1) break;
      case PRIMARY:
       printf("\nThe PRIMARY command\n");
       printf("SYNTAX:\tPRIMARY [residue]+ END\n");
       printf("FUNCTION:\t    \n\n");
       if (ycom != -1) break;
      case READ:
       printf("\nThe READ command\n");
       printf("SYNTAX:\tREAD coord [ref] [format] [<filename>]\n");
       printf("SYNTAX:\tREAD charges [<filename>]\n");
       printf("SYNTAX:\tREAD internal [all|bond|nbond|angle|torsion|rotbond]\n");
       printf("SYNTAX:\t              [<filename>]\n");
       printf("SYNTAX:\tREAD param [<filename>]\n");
       printf("FUNCTION:\t    \n\n");
       if (ycom != -1) break;
      case RECORD:
       printf("\nThe RECORD command\n");
       printf("SYNTAX:\tRECORD [phi=<rotno>] [phib=<value>] [phie=<value>]\n");
       printf("SYNTAX:\t       [psi=<rotno>] [psib=<value>] [psie=<value>]\n");
       printf("SYNTAX:\t       [phiinc=<rotno>] [psiinc=<value>]\n");
       printf("SYNTAX:\t       [recoptc=<value>] [noise|silent]\n");
       printf("SYNTAX:\t       [<filename>]\n");
       printf("FUNCTION:\t    \n\n");
       if (ycom != -1) break;
      case SET:
       printf("\nThe SET command\n");
       printf("SYNTAX:\tSET id <string>\n");
       printf("SYNTAX:\tSET titel \"<string>\"\n");
       printf("FUNCTION:\t    \n\n");
       if (ycom != -1) break;
      case TIME:
       printf("\nThe TIME command\n");
       printf("SYNTAX:\tTIME\n");
       printf("FUNCTION:\t    \n\n");
       if (ycom != -1) break;
      case WRITE:
       printf("\nThe WRITE command\n");
       printf("SYNTAX:\tWRITE internal [all|bond|nbond|angle|torsion|rotbond]\n");
       printf("SYNTAX:\t               [<filename>]\n");
       printf("SYNTAX:\tWRITE param [<filename>]\n");
       printf("FUNCTION:\t    \n\n");
       break;
      otherwise:
       printf("\nThe ????? command\n");
   }

} /* End of phelp */

