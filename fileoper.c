/*
 * This module contains source code for generalised file-operations.
 * Written by Soren Balling Engelsen, INRA-93.
 */
#include "fileoper.h"

/* Generalized routine for opening of files */
FILE *fileopen(char *name, char *mode)
{
   FILE *fp = fopen(name, mode);   /* open the file */

   if (fp == NULL)
   {  fprintf(stderr, "Can't open %s for ", name);
      switch(mode[0])
      {  case 'r': fprintf(stderr, "reading\n");
                   break;
         case 'w': fprintf(stderr, "writing\n");
                   break;
         case 'a': fprintf(stderr, "appending\n");
                   break;
          default: fprintf(stderr, "strange mode\n");
                   break;
      }
   }
   return fp;
}

/* Generalized routine for closing of files */
void
fileclose(char *name, FILE *fp)
{
   if (fp != NULL && fclose(fp) == EOF)
      fprintf(stderr, "Error closing %s\n", name);
}
