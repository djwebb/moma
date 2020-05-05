      subroutine setvbc (ic, jc)
c
c @(#) SCCS module: setvbc.F  version: 1.7
c     Creation date: 06/30/95
c
c==================================================================
c
c     compute vertical boundary conditions for clinic & tracer
c
c==================================================================
c
#include "def_slave.h"
#include "param.h"
c#include "cvbc.h"
#include "grdvar.h"
#include "scalar.h"
#include "slabs.h"
#include "timelv.h"
#include "levind.h"
c
c	TASK COMMON /cvbc/ smf(2),bmf(2),stf(nt),btf(nt)
	COMMON /cvbc/ smf(2),bmf(2),stf(nt),btf(nt)
c-----------------------------------------------------------------
c     get test case surface temperatures, salinities and wind 
c     stresses (they vary with latitude).  test assumes nt=2,
c     (potential temperature and salinity)
c------------------------------------------------------------------
c
      call bcest(stf,smf, ic, jc)
c  convert salinity to model units.
      stf(2) = (stf(2)-35.0)*0.001
c
c------------------------------------------------------------------
c     set no flux condition for tracer at  bottom.
c     relax to "bcest" values with 50 day timescale at top
c------------------------------------------------------------------
c
      gamma = secday/50.0
      do 100 m=1,nt
        stf(m) = gamma*(stf(m)-t(1,ic,jc,m,nm))*dz(1)
        btf(m) = c0
100   continue
c
c------------------------------------------------------------------
c     set bottom stresses.  bmf is stress acting on the ocean floor
c------------------------------------------------------------------
c
      kz = kmu(ic,jc)
      if (kz .ne. 0) then
        uvmag    = sqrt(u(kz,ic,jc,nm)**2 +v(kz,ic,jc,nm)**2)
        bmf(1) = cdbot*u(kz,ic,jc,nm)*uvmag
        bmf(2) = cdbot*v(kz,ic,jc,nm)*uvmag
      else
        bmf(1) = c0
        bmf(2) = c0
      endif
c
      return
      end
