/*
 * This module contain basic I/O rutines for interacting with the 
 * user and with files.
 * Written by Soren Balling Engelsen, INRA-93.
 */
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stddef.h>
#include "defs.h"
#include "sinput.h"

extern BOOLEAN debug;

/* Routine to read a line from stdin or from a file 'fp' */
int
get_line(FILE *fp, char *startptr, int max)
{
   register int c;
   register char *ptr = startptr;
   register char *endptr = startptr + max - 1;

   while ((c = getc(fp)) NOT '\n' AND c != EOF)
      if (ptr < endptr)
         *ptr++ = c;
   *ptr = '\0';
   
   if (debug) printf("GETLINE: %s\n", startptr);

   return (c == EOF) ? -1 : ptr - startptr;
}


/* Routine to write a line to 'stdout' or to a file 'fp' */
void
putline(FILE *fp, char *startptr)
{
   register char *ptr = startptr;

   for (; *ptr != '\0'; ptr++)
      putc(*ptr, fp);
   putc('\n', fp);

}

/* 
 * Routine to place a 'prompt' on screen and get a 'name' (ex.
 * session ID or molecular ID ).
 */
int
enter_name(char *prompt, char *name)
{
   int c;
   char *endptr;     /* points to the last possible character */
   char *startptr;   /* points to the first character */

   printf("%s: ", prompt);
   startptr = name;
   endptr = name + MAXLEN;
   while ((c = getchar()) != '\n' && c != EOF)
      if (name < endptr)
         *name++ = c;
   *name = '\0';
   
   if (debug) printf("ENTER_NAME: [%d] %s\n\n",
                       (name - startptr), startptr);

   return (c == EOF) ? -1 : name - startptr;
}

/*
 * Function that requests a single character input after
 * printing a prompt on the screen.
 */
BOOLEAN
yesorno(char *prompt)
{
   int c;
   
   printf("\n%s: ", prompt);
   while ((c = getchar()) != 'Y' AND c != 'N')
      putchar('\a');
   return(c == 'Y');

}

/* Function that reads a value from the user */
enum error
getvalue(char *prompt, enum itype type, union value *up)
{
   char line[MAXLEN + 1];
   char *cp;
  
   printf("%s", prompt);
   if (get_line(stdin, line, MAXLEN) == -1)
      return NOINPUT;
   switch(type)
   {  case INTEGER: if (sscanf(line, "%d", &(up->i)) != 1)
                       return BADNUMBER;
                    break;
      case REAL:    if (sscanf(line, "%lf", &(up->d)) != 1)
                       return BADNUMBER;
                    break;
      case STRING:  if ( ( up->s = strdup(line) ) == NULL)
                       return NOSPACE;
                    break;
      default:      return BADTYPE;
   }

   return NOERROR;     /* succes */
}

/* Function that print error on stdout */
void
printerror(enum error code)
{
   switch(code)
   {  case NOINPUT:   putline(stdout, "Empty input line.");
                      break;
      case NOSPACE:   putline(stdout, "Out of storage space.");
                      break;
      case BADNUMBER: putline(stdout, "Invalid numeric value.");
                      break;
      case BADTYPE:   putline(stdout, "Internal errormessage; bad type.");
                      break;
      default:        putline(stdout, "Unknown error.");
                      break;
   }
}

/* Function that will prompt for and take an INTEGER from stdin */
int
getinteger(char *prompt)
{
   union value x;
   enum error e;

   while ((e = getvalue(prompt, INTEGER, &x)) != NOERROR)
      printerror(e);
   putchar('\n');

   return x.i;
}
 
/* Function that will prompt for and take an DOUBLE from stdin */
double
getreal(char *prompt)
{
   union value x;
   enum error e;

   while ((e = getvalue(prompt, REAL, &x)) != NOERROR)
      printerror(e);
   putchar('\n');

   return x.d;
}
 

/* Function that will prompt for and take an STRING from stdin */
char
*getstring(char *prompt)
{
   union value x;
   enum error e;

   while ((e = getvalue(prompt, STRING, &x)) != NOERROR)
      printerror(e);
   putchar('\n');

   return x.s;
}

/*
 * This simple routine will print a line of 'width' characters on standard
 * output followed by an empty line.
 *
 * wline('\n',0) shifts one lines
 * wline('\n',n) shifts n+1 lines
 */
void
wline(char ctyp, unsigned short width)
{
   int    i;

   for ( i = 0; i < width; i++) 
      putchar(ctyp);
}
