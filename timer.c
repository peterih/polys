/*
 * This module contains some rutines for obtaining and manipulations
 * with the system clock.
 * Written by Soren Balling Engelsen, INRA-93.
 */
#include <stdio.h>
#include <sys/types.h>
#include <time.h>

/* cput returns elapsed time in seconds since the first call */
double
cput()
{
   double sec;

   sec = 0.000001*(double) clock();
   return(sec);
}

/* */
void
dtime(int *year, int *month, int *day, int *hour, int *min, int *sec)
{
   struct tm *t;
   long clock;

   time(&clock);
   t = localtime(&clock);
   *year = t->tm_year;
   *year += 1900; /* Because tm_year returns year since 1900 */
   *month= t->tm_mon + 1; /* Converts range from 1 to 12 */
   *day  = t->tm_mday;
   *hour = t->tm_hour;
   *min  = t->tm_min;
   *sec  = t->tm_sec;
}

/* */
void
dclock(double *t)
{
   *t = 0.000001*(double) clock();
}

/* */
double
etime(double t)
{
   double ptime;
   
   ptime = 0.000001*(double) clock();

   return(ptime - t); 
}

/* */
void
ptime(double t)
{
   printf("Elapsed CPU time: %1.2f sec.\n", t = etime(0.0));
}

/* */
void
pdate(double t)
{
   int yy, mm, dd, hh, min, ss;

   dtime(&yy,&mm,&dd,&hh,&min,&ss);
   printf("Now:  %d:%02d on %d-%d-%d\n", hh, min, dd, mm, yy);
}
