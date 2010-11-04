/*
 * Written by Soren Balling Engelsen, INRA-93.
 */
#include "polys.h"
#include "extern.h"

void ReadParam(void)
{
   int tstart;

   tstart = cput();
   read_bp(M.nbond);
   read_ap(M.ntheta);
   read_tp(M.nphi);
   read_cp(M.nphi, &M.ncon);
   read_nbp(M.nat, M.nbond);
   tpara += cput() - tstart;
}
