/**********************************************************************
	"sinput.h" Header file for 
SBE93-INRA
**********************************************************************/
#include <stdio.h>
#include "defs.h"

enum itype {INTEGER, REAL, STRING};
enum error {NOERROR, NOINPUT, NOSPACE, BADNUMBER, BADTYPE};

union value          /* can hold a integer, a double or a pointer to a char */
{  int    i;
   double d;
   char   *s;
};

int enter_name(char *prompt, char *name);
int get_line(FILE *fp, char *startptr, int nmax);
void putline(FILE *fp, char *startptr);
BOOLEAN yesorno(char *prompt);
int getinteger(char *prompt);
double getreal(char *prompt);
char *getstring(char *prompt);
void wline(char ctyp, unsigned short width);

/* End of file */
