/*
 * Header file for generalised file operations.
 * SBE93-INRA
 */
#include <stdio.h>
#include <stddef.h>

extern FILE *fileopen(char *name, char *mode);
extern void fileclose(char *name, FILE *fp);
