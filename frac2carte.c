/**********************************************************************
  "frac2carte.c"

  Written by Soren Balling Engelsen, INRA-93.
**********************************************************************/
#include <math.h>
#include "polys.h"
#include "extern.h"

/**********************************************************************
**********************************************************************/
void frac2carte(Crystal UC, BOOLEAN ref)

{  int      i;
   double   ang, vol, cosang[3], sinang[3];
   Vector3  av, avo, avfrac, avcart;

   printf("\n\nConverting fractional coordinates to cartesians\n\n");
   printf("    a = %7.4lf     b = %7.4lf      c = %7.4lf\n\n",
           UC.a, UC.b, UC.c);
   printf("alpha = %7.4lf  beta = %7.4lf  gamma = %7.4lf\n\n",
           UC.alpha, UC.beta, UC.gamma);


   ang = DEGtoRAD*UC.alpha;
   cosang[0] = cos(ang);
   sinang[0] = sin(ang);
   ang = DEGtoRAD*UC.beta;
   cosang[1] = cos(ang);
   sinang[1] = sin(ang);
   ang = DEGtoRAD*UC.gamma;
   cosang[2] = cos(ang);
   sinang[2] = sin(ang);

   /* calculate unit cell volume */
   vol = UC.a*UC.b*UC.c;
   vol *= sqrt(1.0 + 2.0*cosang[0]*cosang[1]*cosang[2] - 
          cosang[0]*cosang[0]-cosang[1]*cosang[1]-cosang[2]*cosang[2]);
   printf("Unit Cell Volume: %15.5lf A^3\n\n", vol);

   CTF[0][0] = UC.a;
   CTF[0][1] = UC.b * cosang[2];
   CTF[0][2] = UC.c * cosang[1];
   CTF[0][3] = 0.0;
   CTF[1][0] = 0.0;
   CTF[1][1] = UC.b * sinang[2];
   CTF[1][2] = UC.c * (cosang[0] - cosang[1]*cosang[2])/sinang[2];
   CTF[1][3] = 0.0;
   CTF[2][0] = 0.0;
   CTF[2][1] = 0.0;
   CTF[2][2] = vol / (UC.a*UC.b*sinang[2]);
   CTF[2][3] = 0.0;
   CTF[3][3] = 1.0;

   printf("Transformation Matrix:\n\n");
   matprint(CTF);

   avfrac.x = 0.0;
   avfrac.y = 0.0;
   avfrac.z = 0.0;
   avcart.x = 0.0;
   avcart.y = 0.0;
   avcart.z = 0.0;
   for (i=0; i<M.nat; i++)
   {  if (ref == FALSE)
        Ar[i] = A[i].pos;
      av = Ar[i];
      avfrac.x += av.x;
      avfrac.y += av.y;
      avfrac.z += av.z;
      TransfV(av, CTF, &avo);
      avcart.x += avo.x;
      avcart.y += avo.y;
      avcart.z += avo.z;
      if (ref == TRUE)
         Ar[i] = avo;
      else
         A[i].pos = avo;
   }

   printf("Geometrical mean values of molecule:\n\n");
   avfrac.x /= M.nat;
   avfrac.y /= M.nat;
   avfrac.z /= M.nat;
   avcart.x /= M.nat;
   avcart.y /= M.nat;
   avcart.z /= M.nat;
   printf("<a> = %7.5lf  <b> = %7.5lf  <c> = %7.5lf\n\n",
           avfrac.x,     avfrac.y,     avfrac.z);
   printf("<x> = %7.5lf  <y> = %7.5lf  <z> = %7.5lf\n\n",
           avcart.x,     avcart.y,     avcart.z);

} /* End of frac2carte */

/* End of file */
