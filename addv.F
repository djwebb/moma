      subroutine addv( ic, jc)
c
c @(#) SCCS module: addv.F  version: 1.10
c     Creation date: 08/09/96
c
c==================================================================
c
c  add barotropic to baroclinic velocities
c
c==================================================================
c
#include "def_slave.h"
#include "param.h"
c
#include "frees.h"
#include "levind.h"
#include "timelv.h"
#include "slabs.h"
c
      ubar=u0(ic,jc,np0)
      vbar=v0(ic,jc,np0)
      do 100 k=1,kmu(ic,jc)
        u(k,ic,jc,np) = u(k,ic,jc,np) + ubar
        v(k,ic,jc,np) = v(k,ic,jc,np) + vbar
 100  continue
c
      return
      end
