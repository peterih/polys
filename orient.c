/*
 * Code for orienting the polysaccharide in space.
 * Written by Soren Balling Engelsen, INRA-93.
 */
#include <stdio.h>
#include "defs.h"
#include "typs.h"
#include "matrix.h"
#include "extern.h"

/*
 * Transform:
 *     atomnumber no  to origo 
 *     atomnumber nx  to the X-axis
 *     atomnumber nxy to the XY-plane 
 * The variables "ia" and "oa" should contain the range of atom
 * indicies to be transformed.
 */
void
Orient(int res1, int res2, int no, int nx, int nxy)
{
   int i;
   Vector3 vo, vx, vxy;
   double  da, db, dc;

   if (debug)
   {  printf("ORIENT: The resiude interval [%d -> %d[\n",res1,res2);
      printf("ORIENT: In origo:    %s[%d]\n", A[no].lab, no);
      printf("ORIENT: On x-axis:   %s[%d]\n", A[nx].lab, nx);
      printf("ORIENT: In xy-plane: %s[%d]\n", A[nxy].lab, nxy);
   }
   /* first translation to origo */
   vo = A[no].pos;
   for (i=0; ((i<M.nat) AND (A[i].res <= res2)); i++)
      if (A[i].res >= res1)
      {  A[i].pos.x -= vo.x;
         A[i].pos.y -= vo.y;
         A[i].pos.z -= vo.z;
      }
   /* then rotate A[nx] to the X-axis */
   vx = A[nx].pos;

   /* and if needed A[nxy] to the XY-plane */
   vxy = A[nxy].pos;

   crossprod(vx, vxy, &vo);
   normalize(&vo);
   normalize(&vx);
   crossprod(vo, vx, &vxy);
   for (i=0; ((i<M.nat) AND (A[i].res <= res2)); i++)
      if (A[i].res >= res1)
      {  da = 0.0; db = 0.0; dc = 0.0;
         da = dotprod(vx, A[i].pos);
         db = dotprod(vxy, A[i].pos);
         dc = dotprod(vo, A[i].pos);
         A[i].pos.x = da;
         A[i].pos.y = db;
         A[i].pos.z = dc;
      }

}
