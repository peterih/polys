/* TOLE = 0.02 used for phimin */

int    ncom=0;    /* defining declarations */
double *pcom=0, *xicom=0, (*nrfunc)();

void linmin(double p[], double xi[], int n, double tole, 
             double *fret, double (*func)())

{  int    j;
   double xx,xmin,fx,fb,fa,bx,ax;
   double brent(),f1dim(),*dvector();
   void   mnbrak(),free_dvector();

   ncom=n;
   pcom=dvector(1,n);
   xicom=dvector(1,n);
   nrfunc=func;
   for (j=1;j<=n;j++) 
   {  pcom[j]=p[j];
      xicom[j]=xi[j];
   }
   ax=0.0;
   xx=1.0;
   bx=2.0;
   mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
   *fret=brent(ax,xx,bx,f1dim,tole,&xmin);
   for (j=1;j<=n;j++) 
   {  xi[j] *= xmin;
      p[j] += xi[j];
   }
   free_dvector(xicom,1,n);
   free_dvector(pcom,1,n);

} /* End of linmin */ 
