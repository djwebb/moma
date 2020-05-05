      subroutine step 
c
c @(#) SCCS module: step.F  version: 1.10
c      Creation date: 08/09/96
c
c=======================================================================
c
c     step is called once per timestep. it includes all the main 
c     loops over ic and jc and calls to the main routines.
c
c=======================================================================
c234567890c234567890c234567890c234567890c234567890c234567890c23456789012
c
#include "def_slave.h"
#include "param.h"
c
#include "scalar.h"
#include "switch.h"
#include "ctmngr.h"
#include "iounit.h"
#include "mesdta.h"
#include "timelv.h"
#include "slabs.h"
#include "frees.h"
c
#include "cdiag.h"
#include "levind.h"
#ifdef inline
# include "chmix.h"
# include "grdvar.h"
# define inline1
#endif
c
c-----------------------------------------------------------------------
c    update pointers for new value of itt.  
c    nnp, nnc and nnm are not changed during a timestep
c    np, nc, nm may be modified during a forward or euler backward 
c    timestep.
c-----------------------------------------------------------------------
c
      nnc = np
      nnm = nc
      nnp = nm
      np  = nnp
      nc  = nnc
      nm  = nnm
c
c-----------------------------------------------------------------------
c       adjust various quantities for normal/mixing timesteps
c-----------------------------------------------------------------------
c
      mxpas2 = .false.
      eots   = .true.
      if (mixts) then
        if (eb) eots = .false.
        nm = nnc
        c2dtts = dtts
        c2dtuv = dtuv
      else
        c2dtts = c2*dtts
        c2dtuv = c2*dtuv
      endif
c
c-----------------------------------------------------------------------
c       return here for second pass of euler backward timestep
c-----------------------------------------------------------------------
c
 100  continue
#ifdef presetp
c
c-----------------------------------------------------------------------
c  
c     precalculate the baroclinic part of the pressure field 
c     for use by subroutine clinic
c
c-----------------------------------------------------------------------
c  
CFPP$ CNCALL
c     do 200 jc=1,jt
c     do 200 ic=1,it
      do 200 l = 1,ncore
        call setp (ia(l),ja(l))
 200  continue 
# ifdef SYNC_3
      do 210 l = 1,nouter
        call setp(iout(l),jout(l))
 210  continue
# endif
#endif
c
c-----------------------------------------------------------------------
c      main baroclinic timestep loop over interior grid cells
c-----------------------------------------------------------------------
c
CFPP$ CNCALL
      call sendposn(14)
      do 300 l = 1, ncore
#ifndef SYNC_3
c
c-----------------------------------------------------------------------
c  check all boundary data has arrived
c-----------------------------------------------------------------------
c
        if(l.eq.ninner+1.and.nproc.gt.1)then
 240      continue
          if(nfalse(1).gt.0)then
            call s_recv(MSG_BOUND1)
            goto 240
          endif
          call initslvd(1)
          call s_send(MSG_CLR_BD1)
# ifdef presetp
          do 250 ll=1,nouter
            call setp(iout(ll),jout(ll))
 250      continue
# endif
        endif
c
c-----------------------------------------------------------------------
c  call baroclinic subroutines
c-----------------------------------------------------------------------
#endif
c
        ic = ia(l)
        jc = ja(l)
        call setvbc (ic,jc)
        call clinic (ic,jc)
        call tracer (ic,jc)
 300  continue
c
c-----------------------------------------------------------------------
c      run free surface model (except during the second part
c      of a baroclinic euler backwards timestep).
c      first initialise pointers 
c-----------------------------------------------------------------------
c
      call sendposn(20)
      if(.not.mxpas2)then
        ntbtp0 = 0
        c2dtbt = dtbt
#ifdef free_eb
        ntbt2 = ntbt
#else
        ntbt2 = 2*ntbt
        do 320 l = 1, ncore
          ic = ia(l)
          jc = ja(l)
          do 320 m=1,3
            freeav(m,ic,jc)=c0
 320    continue
#endif
        do 600 itbt=1,ntbt2
          nnc0 = np0
          nnm0 = nc0
          nnp0 = nm0
          np0  = nnp0
          nc0  = nnc0
#ifdef free_eb
          nm0 = nnc0
#else
          if(itbt.eq.1)then
            nm0 = nnc0
          else
            nm0 = nnm0
          endif
#endif
c
c-----------------------------------------------------------------------
c      On the first barotropic timestep or if "free_eb" is set
c      use euler backward. this requires two passes.
c-----------------------------------------------------------------------
c
          do 590 itbtp=ntbtp0,1
c
c-----------------------------------------------------------------------
c      main free surface model loop to carry out a
c      partial timestep for each model point
c-----------------------------------------------------------------------
c
CFPP$ CNCALL
          do 400 l = 1, ncore
#ifndef SYNC_2
            if(l.eq.ninner+1.and.nproc.gt.1)then
  360         continue
              if(nfalse(2).gt.0)then
                call s_recv(MSG_BOUND2)
                goto 360
              endif
c
c-----------------------------------------------------------------------
c  initialise slave arrays and timestep pointers for 2-D data 
c  then send 2-D ready to receive
c-----------------------------------------------------------------------
c
              call initslvd(2)
              call s_send(MSG_CLR_BD2)
            endif
#endif
#ifdef inline
#         include "frees.F"
#else
          call frees(ia(l),ja(l))
#endif
 400      continue
          if(nproc.gt.1) then
#ifdef SYNC_2
c
c-----------------------------------------------------------------------
c  initialise slave arrays and timestep pointers for 2-D data
c  then send 2-D ready to receive
c-----------------------------------------------------------------------
c
            call initslvd(2)
            call s_send(MSG_CLR_BD2)
#endif
c
c-----------------------------------------------------------------------
c  wait for neighbours permissions then send 2-D boundary data
c-----------------------------------------------------------------------
c
 420        continue
            if(lbtrp)then
              call initslvc(2)
              call s_send(MSG_BOUND2)
            else
              call s_recv(MSG_CLR_BD2)
              goto 420
            endif
#ifdef SYNC_2
c
c-----------------------------------------------------------------------
c    receive 2-d boundary data to outer halo points
c-----------------------------------------------------------------------
c
440         call s_recv(MSG_BOUND2)
            if (nfalse(2) .gt. 0) go to 440
#endif
          endif
c
c-----------------------------------------------------------------
c  set cyclic boundary conditions for the free surface model
c-----------------------------------------------------------------
c
          if(icycl.eq.1) then
            do 500 jc=1,jt
              h0( 1,jc,np0) = h0(itm1,jc,np0)
              h0(it,jc,np0) = h0(   2,jc,np0)
              u0( 1,jc,np0) = u0(itm1,jc,np0)
              u0(it,jc,np0) = u0(   2,jc,np0)
              v0( 1,jc,np0) = v0(itm1,jc,np0)
              v0(it,jc,np0) = v0(   2,jc,np0)
500         continue
          endif
c
c-----------------------------------------------------------------------
c     reset pointers at end of first pass of the free surface model
c-----------------------------------------------------------------------
c
#ifdef free_eb
          if(itbtp.eq.0)then
            nc0 = nnp0
            np0 = nnm0
          else
            nc0 = nnc0
            nm0 = nnp0
          endif
#else
          if(itbt.eq.1)then
            if(itbtp.eq.0)then
              nc0 = nnp0
              np0 = nnm0
            else
              nc0 = nnc0
              nm0 = nnp0
              ntbtp0 = 1
              c2dtbt = c2*dtbt
            endif
          endif
#endif
 590      continue
 600    continue
        itbt = 0
      endif
c
c-----------------------------------------------------------------
c  end of free surface model
c  now add barotropic velocities to baroclinic velocities
c-----------------------------------------------------------------
c
      call sendposn(21)
CFPP$ CNCALL
c     do 700 jc=2,jtm1
c     do 700 ic=2,itm1
      do 700 l = 1, ncore
        call addv(ia(l),ja(l))
 700  continue
c
      if(nproc.gt.1) then
#ifdef SYNC_3
c
c-----------------------------------------------------------------------
c  initialise slave arrays and timestep pointers for 2-D data
c  then send 2-D ready to receive
c-----------------------------------------------------------------------
c
          call initslvd(1)
          call s_send(MSG_CLR_BD1)
#endif
c        
c-----------------------------------------------------------------
c     send 3-d boundary data to other processors
c-----------------------------------------------------------------
c
760       continue
          if(lbcln) then 
            call initslvc(1)
            call s_send(MSG_BOUND1)
          else 
            call s_recv(MSG_CLR_BD1) 
            goto 760 
          endif
#ifdef SYNC_3
c        
c----------------------------------------------------------------
c   receive 3-d data from other processors
c----------------------------------------------------------------
c
780       call s_recv (MSG_BOUND1)
          if (nfalse(1) .gt. 0) goto 780
#endif
      endif
c
c-----------------------------------------------------------------
c  set cyclic boundary conditions for the baroclinic model
c-----------------------------------------------------------------
c
      if(icycl.eq.1) then
        do 800 jc=1,jt
        do 800 k =1,km
          u(k,  1,jc,np) = u(k, itm1,jc,np)
          u(k, it,jc,np) = u(k,    2,jc,np)
          v(k,  1,jc,np) = v(k, itm1,jc,np)
          v(k, it,jc,np) = v(k,    2,jc,np)
          do 800 n=1,nt
            t(k, 1,jc,n,np) = t(k,itm1,jc,n,np)
            t(k,it,jc,n,np) = t(k,   2,jc,n,np)
 800    continue
      endif
c
c-----------------------------------------------------------------
c     if this is the end of the first pass of an euler backward
c     timestep then set the pointers for the second pass.
c-----------------------------------------------------------------
c
      if (mixts.and.eb)then
        eots   = .true.
        nc = nnp
        np = nnm
        mixts  = .false.
        mxpas2 = .true.
        go to 100
      endif
c
c-----------------------------------------------------------------
c  collect timestep statistics
c-----------------------------------------------------------------
c
      call sendposn(32)
      if(prntsi.and.eots) then
        itttsi = itt
        ektot  = c0
        do 900 n = 1,nt
          tddt(n)  = c0
          dtabs(n) = c0
          tvar(n)  = c0
 900    continue
c       do 910 jc=2,jtm1
c       do 910 ic=2,itm1
        do 910 l = 1, ncore
          call diag(ia(l),ja(l))
 910    continue
c
c  set 'ready to send tsi information' flag
c
        ltsiw = .true.
      endif
c
c-----------------------------------------------------------------------
c       if this is the end of either a forward or euler backward
c       timestep then reset the pointers
c-----------------------------------------------------------------------
c
      if(mxpas2)then
        nc = nnc
        nm = nnp
      endif
      if(mixts)then
        nm = nnm
      endif
c
      return
      end
#ifdef inline1
# undef inline1
#endif
