      subroutine setp (ic, jc)
#ifdef presetp
c
c @(#) SCCS module: setp.F  version: 1.1
c     Creation date: 10/13/97
c
c==================================================================
c
c     precalculate the baroclinic part of the pressure field 
c     for use by subroutine clinic
c
c     p(1,...) statement moved within loop 100 to overcome bug
c     with Silicon Graphics at level 3 optimisation
c
c==================================================================
c
#include "def_slave.h"
#include "param.h"
#include "grdvar.h"
#include "scalar.h"
#include "slabs.h"
#include "timelv.h"
#include "levind.h"
#include "mesdta.h"
c
      kmc = kmt(ic,jc)
      if(kmc.gt.0)then        
        call state (t(1,ic,jc,1,nc), t(1,ic,jc,2,nc),  
     &              temp, kmc ,0)
        fxa = p5*grav
        do 100 k=1,kmc
          if(k.eq.1)then
            p(1,ic,jc) = temp(1)*grav*dzw(0)
          else
            p(k,ic,jc) = p(k-1,ic,jc) 
     &                       + dzw(k-1)*fxa*(temp(k)+temp(k-1))
          endif
 100    continue
      endif
c
c------------------------------------------------------------------
c   set cyclic conditions
c------------------------------------------------------------------
c
      if(icycl.eq.1)then
        if(ic.eq.itm1)then
          do 200 k=1,kmc
 200      p(k,1,jc) = p(k,ic,jc)
        endif
        if(ic.eq.2)then
          do 210 k=1,kmc
 210      p(k,it,jc) = p(k,ic,jc)
        endif
      endif
#endif
      return
      end
