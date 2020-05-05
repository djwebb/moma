      program moma
c
c==================================================================
c  Program MOMA,     Version 1.12
c
c  Free surface ocean model for use with array processors.
c
c  Copyright 1993.  D.J. Webb.  
c                   Institute of Oceanographic Sciences,
c                   Godalming, Surrey GU8 5UB, U.K..
c  
c @(#) SCCS information - module:  moma.F     
c @(#)                    version: 1.12    date: 09/22/93
c
c
c  This code is a free surface version of the Bryan-Cox-Semtner
c  code reorganised for efficient running on array processor
c  computers.  The basic finite difference scheme used in the
c  model is described in Bryan (1969), Semtner (1974) and
c  Cox (1984).  The present code follows most closely
c  the format of the GFDL Modular Ocean Model (Pacanowski et al,
c  1990), and makes use of many of the subroutines and include
c  files of that model.  More detailed information about the 
c  present code is given in Webb (1993).
c
c  The code may be freely adapted and used without charge.
c
c  For further information on the model contact the author.
c  e-mail:   djw@unixa.nerc-wormley.ac.uk
c  telemail: ios.wormley
c
c  
c  The main differences from the modular ocean model code code are:
c
c  1.  Collection of all 'array processor' loops (i.e. loops over the
c      horizontal indices ic and jc) in subroutine step.  When an array
c      processor is used, code should be added to this subroutine to 
c      partition the horizontal index (ic, jc) ranges between the
c      different processors so that they each have a similar workload.
c  2.  Revised common block structures.  All variables are stored
c      in core.
c  3.  Removal of all 'slab' optimisation code in the program and
c      the introduction of code designed to optimise inner loop 
c      calculations using the vertical index k.
c  4.  To simplify the development and testing of the array processor
c      code the following features of the moma code were removed:
c        (a) Diagnostic calculations in routines clinic and tracer.
c        (b) Moma code options.
c      It should now be straightforward (but possibly time consuming)
c      to add the features back into the present code if they are 
c      required:
c  5.  The stream function code is removed and replaced by a free-
c      surface code.  The latter is similar to Killworth, Stainforth
c      Webb and Patterson (1989).  For efficiency, the mean horizontal
c      velocity is now used instead of horizontal transport and
c      the viscous terms are calculated in routine clinic.
c  6.  In the baroclinic momentum equation, a revised horizontal
c      advection scheme option is included.  The old scheme can still be 
c      used by specifying the flag 'oldadv'
c  7.  Near topography, in the baroclinic momentum equation, a revised
c      vertical advection scheme option is included.
c  8.  An option to precalculate the baroclinic part of the 
c      pressure filed is included using flag 'presetp'.
c
c
c  cpp precompiler options:
c
c     'oldadv'   - use the origonal scheme for the horizontal
c                  advection velocity at velocity points.
c     'presetp'  - precalculate the baroclinic pressure field
c                  before calling clinic
c     'hcomments'- include comments from '*.h' files
c
c  References:
c
c  Bryan, K., 1969:  A numerical method for the circulation of the
c  World Ocean.  Journal of Computational Physics, 4, 347-  .
c
c  Semtner, A.J., 1974:  A general circulation model for the 
c  World Ocean.  UCLA Department of Meteorology Technical Report
c  No. 8, 99pp.
c
c  Cox, M.D., 1984:  A primitive equation, 3-dimensional model of
c  the ocean.  GFDL Ocean Technical Report No.1, Geophysical
c  Fluid Dynamics Laboratory/NOAA, Princeton University, 
c  Princeton N.J., U.S.A..
c
c  Killworth, P.D., Stainforth, D., Webb, D.J. and Paterson, P.M.,
c  1989:  A free surface Bryan-Cox-Semtner model.  Report No. 270,  
c  Institute of Oceanographic Sciences, Wormley, U.K..
c
c  Pacanowski, R.C., Dixon, K., Rosati, A., 1990:  The GFDL Modular
c  Ocean Model 1.0.  Geophysical Fluid Dynamics Laboratory/NOAA, 
c  Princeton University, Princeton, N.J., U.S.A..  (Unpublished 
c  manuscript).
c
c  Webb, D.J., 1993:  An ocean model code for array processor
c  computers.  Internal Document No.324, Institute of Oceanographic 
c  Sciences, Wormley, U.K..
c
c==================================================================
c
#include "param.h"

#include "scalar.h"
#include "switch.h"
#include "coord.h"
#include "ctmngr.h"
#include "grdvar.h"
#include "levind.h"
#include "timelv.h"
#include "slabs.h"

#include "cdiag.h"
#include "chmix.h"
#include "cvbc.h"
#include "cvmix.h"
#include "frees.h"
#include "iounit.h"
#include "versno.h"
c==================================================================
c
      dimension utime(3)
c
      namelist /eddy/   am, ah, fkpm, fkph, cdbot
      namelist /tsteps/ dtts, dtuv, dtbt
      namelist /contrl/ init, fnrest, days, restrt, nmix, eb, ncon, 
     &                  tsi,  dgnstc, snaps, acor
c
c  start timing
c 
      call dtime(utime)
      utime(3)=0.0
c
c     constants usually initialised in block data
c
      model(1) =' Stripped down MOM code for array processor.'
      model(2) =' Module: moma.F.  Version: 1.12. '//
     &          ' Date: 09/22/93'
      stlon = -4.0
      stlat = -72.0
      dxdeg = 4.0
      dydeg = 4.0
c
      dz( 1) = 30.00e2
      dz( 2) = 46.15e2
      dz( 3) = 68.93e2
      dz( 4) = 99.93e2
      dz( 5) = 140.63e2
      dz( 6) = 192.11e2
      dz( 7) = 254.76e2
      dz( 8) = 327.95e2
      dz( 9) = 409.81e2
      dz(10) = 497.11e2
      dz(11) = 585.36e2
      dz(12) = 669.09e2
      dz(13) = 742.41e2
      dz(14) = 799.65e2
      dz(15) = 836.10e2
c
      init = .true.
      first= .true.
      restrt = .false.
      fnrest = ' '
      iorest = 21
      tsi   = 2.0
      dgnstc= 2.0
      snaps = 2.0
c
      dtts = 10800d0
      dtuv = 10800d0
      dtbt = 100d0
c
      am    = 1.0e9 
      ah    = 2.0e7
      fkph  = 20.0
      fkpm  = 1.0
      cdbot = 0.0
      eb    = .true.
      acor  = 0.6
      nmix  = 16
      ncon  = 1
c
c==================================================================
c     begin introductory section
c     open input file and write MOM version information
c==================================================================
c
      open (stdin,file='ocean.in')
      write (stdout,'(/2(25x,a80/)/)') model
c
c------------------------------------------------------------------
c     initialize various quantities and
c     read in run parameters
c------------------------------------------------------------------
c
      pi     = c4*atan(c1)
      radian = c360/(c2*pi)
      omega  = pi/43082.0
      grav =980.6
      radius =6370.e5
c
      read  (stdin, contrl)
      write (stdout,contrl)
      read  (stdin, eddy)
      write (stdout,eddy)
      read  (stdin, tsteps)
      write (stdout,tsteps)
      close (unit=stdin,status='keep')
c
c------------------------------------------------------------------
c     calculate free surface model (barotropic) timestep
c------------------------------------------------------------------
c
      ntbt = max(1.0,dtuv/dtbt)
      if(dtbt.ne.dtuv/ntbt)then
        dtbt = dtuv/ntbt
        write(stdout,'(a)') ' dtuv not integer multiple of dtbt.'
        write(stdout,'(a,i5)') ' dtbt now set to set to:',dtuv
      endif
c
c------------------------------------------------------------------
c     set up model grids in x, y, & z
c------------------------------------------------------------------
c
      call grids
c
c==================================================================
c     end introductory section
c==================================================================
c
c------------------------------------------------------------------
c     set up initial pointers to timestep storage.  nc is current 
c     timestep, np is next and nm is previous timestep.  nm0, nc0
c     and np0 are the corresponding pointers for the free surface
c     model.
c------------------------------------------------------------------
c
      nm = 1
      nc = 2
      np = 3
      nm0 = 1
      nc0 = 2
      np0 = 3
c
c------------------------------------------------------------------
c     is this a start from initial conditions or a restart?
c     in either case, set the number of vertical levels on the "t" 
c     grid
c------------------------------------------------------------------
c
      if (init) then
        call ocn1st
      else
        call rdrest
      endif
c
c------------------------------------------------------------------
c     compute number of vertical levels on the "u" grid
c     and set cyclic boundary conditions
c------------------------------------------------------------------
c
      do 900 i=1,imt
        kmu(i,jmt) = 0
900   continue
      do 1000 j=1,jmtm1
      do 990 i=1,imtm1
        kmu(i,j) =min (kmt(i,j),kmt(i+1,j),kmt(i,j+1),kmt(i+1,j+1))
990   continue
1000  continue
      do 1100 j=1,jmt
        kmu(imt,j) = kmu(2,j)
1100  continue
c
c------------------------------------------------------------------
c     compute area and volume of ocean ("t,s" grid boxes)
c     and set cyclic boundary conditions
c------------------------------------------------------------------
c
      area   = c0
      volume = c0
      ocnp   = 0
c
      do 700 j=2,jmtm1
      do 690 i=2,imtm1
      if (kmt(i,j) .gt. 0) then
        area   = area + cst(j)*dx*dy
        volume = volume + cst(j)*dx*dy*zw(kmt(i,j))
        ocnp   = ocnp + float(kmt(i,j))
      endif
690   continue
700   continue
      write (stdout,9341) area, volume
c
c------------------------------------------------------------------
c       print map of "kmt" levels. 
c------------------------------------------------------------------
c
c
      write (stdout,9402)
      do 1300 ibk=1,imt,40
      isp  = ibk
      iept = ibk + 40 - 1
      if(iept.gt.imt) iept=imt
      write (stdout,'(/, 4x, 40i3)') (ii, ii=isp,iept)
      do 1290 jrev=1,jmt
      j=jmt-jrev+1
      write (stdout,'(1x,i3, 40i3)')j,(kmt(i,j),i=isp,iept)
1290  continue
1300  continue
c
c------------------------------------------------------------------
c     compute depths and reciprocal depths
c------------------------------------------------------------------
c
      do 1400 j=1,jmt
      do 1390 i=1,imt
        h(i,j)  = c0
        hr(i,j) = c0
        if (kmu(i,j) .ne. 0) then
          h (i,j) = zw(kmu(i,j))
          hr(i,j) = c1/zw(kmu(i,j))
        endif
1390  continue
1400  continue
c
c------------------------------------------------------------------
c     initialize various things
c------------------------------------------------------------------
c
      do 2600 j=1,jmt
      do 2590 i=1,imt
      zu(i,j)  = c0
      zv(i,j)  = c0
2590  continue
2600  continue
c
c     initialize the coriolis factor
c
      do 2660 j=1,jmt
      fcor(j) = c2*omega*sine(j)
2660  continue
c
c     initialize the diffusion factors
c
      do 2670 j=2,jmtm1
      bbt(j) = ah*dxr*cstr(j)*dy
      cct(j) = ah*dyr*dx*csu(j  )
      ddt(j) = ah*dyr*dx*csu(j-1)
      bbu(j) = am*dxr*csur(j)*dy
      ccu(j) = am*dyr*dx*cst(j+1)
      ddu(j) = am*dyr*dx*cst(j  )
      ggu(j) = am*(c1-tng(j)*tng(j))/(radius*radius)
      hhu(j) = am*c2*sine(j)/(radius*csu(j)*csu(j))
2670  continue
c
c------------------------------------------------------------------
c     start the time step loop
c------------------------------------------------------------------
c
      do 3400 loop=1,9999999
      first= loop .eq. 1
c
c------------------------------------------------------------------
c       update timestep, set time dependent logical switches
c       to determine program flow for timestep itt, itt-1 and itt-2.
c       note: timestep itt and time refer to the time of the new
c       fields being calculated!
c------------------------------------------------------------------
c
      call tmngr (dtts)
      call step 
c
c-----------------------------------------------------------------
c      print time step information, print diagnostics and archive 
c      as required
c-----------------------------------------------------------------
c
      call dtime(utime)
      utime(3)=utime(1)+utime(2)+utime(3)
      if (prntsi) write (stdout,9602) itt, stamp, ektot,
     &                      dtabs(1), dtabs(2),utime(1),utime(3)
      if (diagts) call pdiag
      if (snapts) call arch
c
c------------------------------------------------------------------
c       end timestepping loop
c------------------------------------------------------------------
c
      if (last) go to 3401
3400  continue
3401  continue
c
c------------------------------------------------------------------
c     save restart at end of run if restrt is true
c------------------------------------------------------------------
c
      if (restrt) call wrrest
c
c------------------------------------------------------------------
c     close all units  (stderr is not closed because it is set equal
c     to stdout in pconst.h).
c------------------------------------------------------------------
c
      write (stdout,9603) 
      close (unit=stdout,status='keep')
      stop
c
9341  format (//,'  Regional & Global ocean statistics:'
     &,/,'  the total ocean surface area (t grid) =',1pe15.8,'cm**2'
     &,/,'  the total ocean volume (t grid)       =',1pe15.8,'cm**3')
9402  format(/t50,'number of levels on "t,s" grid')
9451  format (/' ==== start and end indices for',a17,'====')
9461  format (' j=',i3,5x,5(2i5,10x))
9499  format (/' error => lseg too small for',a15,' indices'
     &        /'          j =',i5,'  lseg + 1 =',i8)
9602  format (1x,'ts=',i7, 1x, a32, ', ke=', 1pe13.6,
     &        ' dtemp=',1pe13.6,' dsalt=',1pe13.6,' times=',0p2f10.2)
9603  format(1x,' ==> END of model run.',//)
      end


      subroutine step 
c
c==================================================================
c
c     step is called once per timestep. it includes all the main 
c     loops over ic and jc and calls to the main routines.
c
c==================================================================
c
c
#include "param.h"
c
#include "scalar.h"
#include "switch.h"
#include "timelv.h"
#include "slabs.h"
#include "frees.h"
c
#include "cdiag.h"
c
c------------------------------------------------------------------
c    update pointers for new value of itt.  
c    nnp, nnc and nnm are not changed during a timestep
c    np, nc, nm may be modified during a forward or euler backward 
c    timestep.
c------------------------------------------------------------------
c
      nnc = np
      nnm = nc
      nnp = nm
      np  = nnp
      nc  = nnc
      nm  = nnm
c
c------------------------------------------------------------------
c       adjust various quantities for normal/mixing timesteps
c------------------------------------------------------------------
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
c------------------------------------------------------------------
c       return here for second pass of euler backward timestep
c------------------------------------------------------------------
c
 100  continue
#ifdef presetp
c
c-----------------------------------------------------------------
c  
c     precalculate the baroclinic part of the pressure field 
c     for use by subroutine clinic
c
c-----------------------------------------------------------------
c  
      do 150 jc=1,jmt
      do 150 ic=1,imt
        call setp (ic,jc)
 150  continue 
#endif
c
c-----------------------------------------------------------------
c      main baroclinic timestep loop over grid cells, from south 
c      to north and from west to east
c      1. set vertical boundary conditions (surface & bottom)
c      2. calculate internal mode velocities
c      3. calculate tracers
c-----------------------------------------------------------------
c
      do 200 jc=2,jmtm1
      do 200 ic=2,imtm1
        call setvbc (ic,jc)
        call clinic (ic,jc)
        call tracer (ic,jc)
 200  continue
c
c-----------------------------------------------------------------
c      run free surface model (except during the second part
c      of a baroclinic euler backwards timestep).
c      first initialise pointers 
c-----------------------------------------------------------------
c
      if(.not.mxpas2)then
        do 600 lt=1,ntbt
          nnc0 = np0
          nnm0 = nc0
          nnp0 = nm0
          np0  = nnp0
          nc0  = nnc0
          nm0  = nnc0
          frpas1=.true.
c
c-----------------------------------------------------------------
c      use an euler backward scheme. this requires two passes.
c-----------------------------------------------------------------
c
          do 500 lb=1,2
c
c-----------------------------------------------------------------
c      main free surface model loop to carry out a
c      partial timestep for each model point
c-----------------------------------------------------------------
c
          do 300 jc=2,jmtm1
          do 300 ic=2,imtm1
            call frees(ic,jc)
 300      continue
c
c-----------------------------------------------------------------
c  set cyclic boundary conditions for the free surface model
c-----------------------------------------------------------------
c
          do 400 jc=2,jmt-1
            h0(  1,jc,np0) = h0(imtm1,jc,np0)
            h0(imt,jc,np0) = h0(    2,jc,np0)
            u0(  1,jc,np0) = u0(imum1,jc,np0)
            u0(imu,jc,np0) = u0(    2,jc,np0)
            v0(  1,jc,np0) = v0(imum1,jc,np0)
            v0(imu,jc,np0) = v0(    2,jc,np0)
400       continue
c
c-----------------------------------------------------------------
c     reset pointers at end of first pass of the free surface model
c------------------------------------------------------------------
c
          if(frpas1)then
            frpas1 = .false.
            nc0 = nnp0
            np0 = nnm0
          endif
 500      continue
 600    continue
      endif
c
c-----------------------------------------------------------------
c  end of free surface model
c  now add barotropic velocities to baroclinic velocities
c-----------------------------------------------------------------
c
      do 700 jc=2,jmtm1
      do 700 ic=2,imtm1
        call addv(ic,jc)
 700  continue
c
c-----------------------------------------------------------------
c  set cyclic boundary conditions for the baroclinic model
c-----------------------------------------------------------------
c
      do 800 jc=2,jmt-1
      do 800 k =1,km
        u(k,  1,jc,np) = u(k,imum1,jc,np)
        u(k,imu,jc,np) = u(k,    2,jc,np)
        v(k,  1,jc,np) = v(k,imum1,jc,np)
        v(k,imu,jc,np) = v(k,    2,jc,np)
        do 800 n=1,nt
          t(k,  1,jc,n,np) = t(k,imtm1,jc,n,np)
          t(k,imt,jc,n,np) = t(k,    2,jc,n,np)
 800  continue
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
      if(prntsi.and.eots) then
        ektot = c0
        do 900 n=1,nt
        tddt(n)  = c0
        dtabs(n) = c0
        tvar(n)  = c0
 900    continue
        do 910 jc=2,jmtm1
        do 910 ic=2,imtm1
        call diag(ic,jc)
 910    continue
        ektot = ektot/volume
        do 920 n=1,nt
        tddt(n)  = tddt(n)/volume
        dtabs(n) = dtabs(n)/volume
        tvar(n)  = tvar(n)/volume
 920    continue
      endif
c
c------------------------------------------------------------------
c       if this is the end of either a forward or euler backward
c       timestep then reset the pointers
c------------------------------------------------------------------
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
 

      subroutine clinic(ic,jc)
c
c==================================================================
c
c     clinic computes internal mode velocities at time level t+1
c     for a vertical column of k cells with i and j values equal to
c     ic and jc.
c
c==================================================================
c
c
#include "param.h"
c
#include "scalar.h"
#include "switch.h"
#include "coord.h"
#include "grdvar.h"
#include "levind.h"
#include "timelv.h"
#include "slabs.h"
#include "frees.h"
c
#include "chmix.h"
#include "cvbc.h"
#include "cvmix.h"
c
      kmc=kmu(ic,jc)
      if(kmc.eq.0)return
      jp = jc+1
      jm = jc-1
      ip = ic+1
      im = ic-1
c
c  find max of all the neighbouring points.
c
#ifdef oldadv
      kmd  = max(kmu(ip,jc),kmu(im,jc),kmu(ic,jp),kmu(ic,jm),kmc)
#else
      kmd  = max(kmu(im,jm),kmu(ic,jm),kmu(ip,jm),
     &           kmu(im,jc),kmc       ,kmu(ip,jc),
     &           kmu(im,jp),kmu(ic,jp),kmu(ip,jp))
#endif
      boxa = dx*csu(jc)*dy
      boxar = c1/boxa
c
c=============================================================
c  calculate horizontal and vertical advective fluxes
c=============================================================
c
      dxp5 = p5*dx
      dyp5 = p5*dy
      do 100 k=1,kmd
#ifdef oldadv
c
c  original horizontall advection scheme
c
      fvn(k) = (v(k,ic,jc,nc)+v(k,ic,jp,nc))*p5*dx*cst(jp)*dz(k)
      fvs(k) = (v(k,ic,jc,nc)+v(k,ic,jm,nc))*p5*dx*cst(jc)*dz(k)
      fue(k) = (u(k,ic,jc,nc)+u(k,ip,jc,nc))*p5*dy*dz(k)
      fuw(k) = (u(k,ic,jc,nc)+u(k,im,jc,nc))*p5*dy*dz(k)
#else
c
c  new horizontal advection scheme
c
      fvjc   = csu(jc)*
     &         (v(k,im,jc,nc)+2.0*v(k,ic,jc,nc)+v(k,ip,jc,nc))
      fuic   = (u(k,ic,jp,nc)+2.0*u(k,ic,jc,nc)+u(k,ic,jm,nc))
c
      fvn(k) = dz(k)*dx*p125*( fvjc + csu(jp)*
     &         (v(k,im,jp,nc)+2.0*v(k,ic,jp,nc)+v(k,ip,jp,nc)))
      fvs(k) = dz(k)*dx*p125*( fvjc + csu(jm)*
     &         (v(k,im,jm,nc)+2.0*v(k,ic,jm,nc)+v(k,ip,jm,nc)))
c
      fue(k) = dz(k)*dy*p125*( fuic +
     &         (u(k,ip,jp,nc)+2.0*u(k,ip,jc,nc)+u(k,ip,jm,nc)))
      fuw(k) = dz(k)*dy*p125*( fuic +
     &         (u(k,im,jp,nc)+2.0*u(k,im,jc,nc)+u(k,im,jm,nc)))
#endif
c
c  new vertical boundary conditions on w at u points
c
      fw(k-1)  = (fuw(k)-fue(k)+fvs(k)-fvn(k))
 100  continue
c
c=============================================================
c    integrate up from the bottom level
c=============================================================
c
      fw(kmd)=c0
      do 210 k=kmd,1,-1
      fw(k-1)=fw(k-1)+fw(k)
 210  continue
#ifdef presetp
c
c-----------------------------------------------------------------
c         compute hydrostatic pressure gradient (baroclinic part)
c-----------------------------------------------------------------
c
      fxa = dx2r*csur(jc)
      fxb = dy2r
      do 300 k=1,kmc
        temp1  = p(k,ip,jp) - p(k,ic,jc)
        temp2  = p(k,ic,jp) - p(k,ip,jc)
        dpdx(k) = (temp1-temp2)*fxa
        dpdy(k) = (temp1+temp2)*fxb
 300  continue
#else
c
c-----------------------------------------------------------------
c     compute density on surrounding 't' points
c-----------------------------------------------------------------
c
      call state (t(1,ic,jc,1,nc), t(1,ic,jc,2,nc), rhoo, 
     &            kmt(ic,jc) ,0)
      call state (t(1,ip,jc,1,nc), t(1,ip,jc,2,nc), rhpo, 
     &            kmt(ip,jc) ,0)
      call state (t(1,ip,jp,1,nc), t(1,ip,jp,2,nc), rhpp, 
     &            kmt(ip,jp) ,0)
      call state (t(1,ic,jp,1,nc), t(1,ic,jp,2,nc), rhop, 
     &            kmt(ic,jp) ,0)
c
c-----------------------------------------------------------------
c         compute hydrostatic pressure gradient (baroclinic part)
c-----------------------------------------------------------------
c
c    1. compute it at the first level
c
      fxa = grav*dzw(0)*dx2r*csur(jc)
      fxb = grav*dzw(0)*dy2r
      temp1  = rhpp(1) - rhoo(1)
      temp2  = rhop(1) - rhpo(1)
      dpdx(1) = (temp1-temp2)*fxa
      dpdy(1) = (temp1+temp2)*fxb
c
c    2. compute the change in pressure gradient between levels
c
      fxa = grav*dx4r*csur(jc)
      fxb = grav*dy4r
      do 300 k=2,kmc
        temp1 = rhpp(k)+rhpp(k-1)- rhoo(k)-rhoo(k-1)
        temp2 = rhop(k)+rhop(k-1)- rhpo(k)-rhpo(k-1)
        dpdx(k) = (temp1-temp2)*fxa*dzw(k-1)
        dpdy(k) = (temp1+temp2)*fxb*dzw(k-1)
 300  continue
c
c=============================================================
c         integrate downward from the first level
c=============================================================
c
      do 400 k=2,kmc
        dpdx(k) = dpdx(k-1) + dpdx(k)
        dpdy(k) = dpdy(k-1) + dpdy(k)
 400  continue
#endif
c
c-----------------------------------------------------------------
c         calculate quantities for the computation of
c         vertical diffusion and advection of momentum
c         - positive for upwards flux of momentum
c-----------------------------------------------------------------
c
      do 500 k=1,kmc-1
        vmf(k,1) = fkpm*(u(k+1,ic,jc,nm) - u(k  ,ic,jc,nm))
     &                     *dzwr(k)*dx*csu(jc)*dy
        vmf(k,2) = fkpm*(v(k+1,ic,jc,nm) - v(k  ,ic,jc,nm))
     &                     *dzwr(k)*dx*csu(jc)*dy
        fwb1(k) = fw(k)*(u(k,ic,jc,nc) + u(k+1,ic,jc,nc))*p5
        fwb2(k) = fw(k)*(v(k,ic,jc,nc) + v(k+1,ic,jc,nc))*p5
 500  continue
c
c-----------------------------------------------------------------
c     set k=0 vertical diffusion to reflect surface wind stress
c     set level "kmc" for bottom drag condition.
c     set k=kmc elements of vertical advection array to zero
c     vertical momentum transfer.  Set k=0 to represent advection
c     loss to region above the free surface.
c     smf is the surface stress acting on the ocean
c     bmf is the bottom stress acting on the ocean floor
c-----------------------------------------------------------------
c
      vmf(0,1)  = -smf(1)*boxa
      vmf(0,2)  = -smf(2)*boxa
      vmf(kmc,1)= -bmf(1)*boxa
      vmf(kmc,2)= -bmf(2)*boxa
      fwb1(0)   = fw(0)*u(1,ic,jc,nc)
      fwb2(0)   = fw(0)*v(1,ic,jc,nc)
      fwb1(kmc) = c0
      fwb2(kmc) = c0
c
c------------------------------------------------------------------
c       set up index lcor for the coriolis term:
c       use tau timestep for explicit trtmnt ... or
c       use tau-1 timestep for implicit treatment with remainder of
c       term to be added later
c------------------------------------------------------------------
c
          if (acor .eq. c0) then
            lcor = nc
          else
            lcor = nm
          endif
c
c------------------------------------------------------------------
c       initialise vertical mean forcing array elements:
c------------------------------------------------------------------
c
      zu(ic,jc) = c0
      zv(ic,jc) = c0
c
c=============================================================
c    main calculation to timestep baroclinic velocity
c=============================================================
c
      do 600 k=1,kmc
c
c u-momentum equation.  separate coriolis term and the rest.
c
      ucor  = fcor(jc)*v(k,ic,jc,lcor)
      urest = -dpdx(k)
c  u-advection
     & +(fuw(k)*(u(k,im,jc,nc) + u(k,ic,jc,nc))*p5
     & - fue(k)*(u(k,ip,jc,nc) + u(k,ic,jc,nc))*p5
c  v-advection
     & + fvs(k)*(u(k,ic,jm,nc) + u(k,ic,jc,nc))*p5
     & - fvn(k)*(u(k,ic,jp,nc) + u(k,ic,jc,nc))*p5
c  w-advection
     & + fwb1(k  ) - fwb1(k-1)
c  x diffusion  
     & +(bbu(jc)*(u(k,ip,jc,nm)-u(k,ic,jc,nm))
     & + bbu(jc)*(u(k,im,jc,nm)-u(k,ic,jc,nm))
c  y diffusion
     & + ccu(jc)*(u(k,ic,jp,nm)-u(k,ic,jc,nm))
     & + ddu(jc)*(u(k,ic,jm,nm)-u(k,ic,jc,nm)) )
     & *dz(k)
c  z diffusion  
     &  + vmf(k  ,1) - vmf(k-1,1) )*boxar*dzr(k)
c
c v-momentum equation.  separate coriolis term and the rest
c
      vcor  = -fcor(jc)*u(k,ic,jc,lcor)
      vrest = -dpdy(k)
c  u-advection
     & +(fuw(k)*(v(k,im,jc,nc) + v(k,ic,jc,nc))*p5
     & - fue(k)*(v(k,ip,jc,nc) + v(k,ic,jc,nc))*p5
c  v-advection
     & + fvs(k)*(v(k,ic,jm,nc) + v(k,ic,jc,nc))*p5
     & - fvn(k)*(v(k,ic,jp,nc) + v(k,ic,jc,nc))*p5
c  w-advection
     & + fwb2(k  ) - fwb2(k-1)
c  x diffusion  
     & +(bbu(jc)*(v(k,ip,jc,nm)-v(k,ic,jc,nm))
     & + bbu(jc)*(v(k,im,jc,nm)-v(k,ic,jc,nm))
c  y diffusion
     & + ccu(jc)*(v(k,ic,jp,nm)-v(k,ic,jc,nm))
     & + ddu(jc)*(v(k,ic,jm,nm)-v(k,ic,jc,nm)) )
     & *dz(k)
c  z diffusion
     & + vmf(k  ,2) - vmf(k-1,2) )*boxar*dzr(k)
c
c***************************************************************
c  the corrections to the advection and diffusion terms that
c  allow for the curvature of the grid near the poles have been
c  left out of this version of the model
c***************************************************************
c
c-----------------------------------------------------------------
c     calculate vertically averaged forcing (without coriolis)
c     and total force at each level (with coriolis)
c-----------------------------------------------------------------
c
      fx = dz(k)
      zu(ic,jc) = zu(ic,jc) + urest*fx
      zv(ic,jc) = zv(ic,jc) + vrest*fx
      u(k,ic,jc,np) = ucor + urest
      v(k,ic,jc,np) = vcor + vrest
 600  continue
c
c-----------------------------------------------------------------
c     calculate vertical mean force 
c-----------------------------------------------------------------
c
      zu(ic,jc) = zu(ic,jc)*hr(ic,jc)
      zv(ic,jc) = zv(ic,jc)*hr(ic,jc)
c
c-----------------------------------------------------------------
c     timestep baroclinic terms
c-----------------------------------------------------------------
c
      if(acor.eq.c0)then
        do 700 k=1,kmc
          u(k,ic,jc,np) = u(k,ic,jc,nm) + c2dtuv*u(k,ic,jc,np) 
          v(k,ic,jc,np) = v(k,ic,jc,nm) + c2dtuv*v(k,ic,jc,np) 
 700    continue
      else
        fx    = c2dtuv*acor*c2*omega*sine(jc)
        detmr = c1/(c1 + fx*fx)
        do 750 k=1,kmc
          t1 = (u(k,ic,jc,np)+fx*v(k,ic,jc,np))*detmr
          t2 = (v(k,ic,jc,np)-fx*u(k,ic,jc,np))*detmr
          u(k,ic,jc,np) = u(k,ic,jc,nm) + c2dtuv*t1 
          v(k,ic,jc,np) = v(k,ic,jc,nm) + c2dtuv*t2
 750    continue
      endif
c
c-----------------------------------------------------------------
c     remove vertical mean velocity
c-----------------------------------------------------------------
c
      sumu = c0
      sumv = c0
      zwr = c1/zw(kmc)
      do 800 k=1,kmc
      sumu = sumu + u(k,ic,jc,np)*dz(k)
      sumv = sumv + v(k,ic,jc,np)*dz(k)
 800  continue
      sumu = sumu*zwr
      sumv = sumv*zwr
      do 900 k=1,kmc
      u(k,ic,jc,np) = u(k,ic,jc,np) - sumu
      v(k,ic,jc,np) = v(k,ic,jc,np) - sumv
 900  continue
      return
      end
 

      subroutine tracer(ic,jc)
c
c=============================================================
c
c     tracer computes tracer quantities for time level t+1
c
c=============================================================
c
#include "param.h"
c
#include "scalar.h"
#include "grdvar.h"
#include "levind.h"
#include "timelv.h"
#include "slabs.h"
c
#include "chmix.h"
#include "cvbc.h"
#include "cvmix.h"
c
c=============================================================
c    outer loops over i and j indices
c=============================================================
c
      kmc= kmt(ic,jc)
      if(kmc.eq.0)return
      jp = jc+1
      jm = jc-1
      ip = ic+1
      im = ic-1
      boxa = dx*cst(jc)*dy
      boxar = c1/boxa
c
c=============================================================
c    set up masks, equal to 0 on sub-seabottom points.
c=============================================================
c
      do 110 k=1,kmc
      maskop(k)=0
      if(k.le.kmt(ic,jp))maskop(k)=1
      maskom(k)=0
      if(k.le.kmt(ic,jm))maskom(k)=1
      maskpo(k)=0
      if(k.le.kmt(ip,jc))maskpo(k)=1
      maskmo(k)=0
      if(k.le.kmt(im,jc))maskmo(k)=1
 110  continue 
c
c=============================================================
c  calculate horizontal and vertical advective fluxes
c=============================================================
c
      dxp5 = p5*dx
      dyp5 = p5*dy
      do 160 k=1,kmc
      fvn(k) = (v(k,ic,jc,nc)+v(k,im,jc,nc))*dxp5*csu(jc)*dz(k)
      fvs(k) = (v(k,ic,jm,nc)+v(k,im,jm,nc))*dxp5*csu(jm)*dz(k)
      fue(k) = (u(k,ic,jc,nc)+u(k,ic,jm,nc))*dyp5*dz(k)
      fuw(k) = (u(k,im,jc,nc)+u(k,im,jm,nc))*dyp5*dz(k)
      fw(k-1)  = (fuw(k)-fue(k)+fvs(k)-fvn(k))
 160  continue

c
c=============================================================
c    integrate up from the lowest layer 
c    fw is positive upwards
c=============================================================
c
      fw(kmc)=c0
      do 180 k=kmc-1,0,-1
      fw(k)=fw(k)+fw(k+1)
 180  continue
c
c=============================================================
c    main calculation loop to timestep tracers
c=============================================================
c
      do 900 m=1,nt
c
c=============================================================
c  calculate vertical fluxes - positive upwards
c=============================================================
c
      do 930 k=1,kmc-1
        vtf(k) = fkph*(t(k+1,ic,jc,m,nm)-t(k  ,ic,jc,m,nm))
     &                 *dzwr(k)*boxa
        fwb1(k) = fw(k)*(t(k,ic,jc,m,nc)+t(k+1,ic,jc,m,nc))*p5
 930  continue
c
c=============================================================
c  surface and bottom advective and diffusive boundary conditions
c  stf is downwards flux into the ocean
c=============================================================
c
      vtf(0)   = -stf(m)*boxa
      vtf(kmc) = c0
      fwb1(0)  = fw(0)*t(1,ic,jc,m,nc)
      fwb1(kmc)= c0
c
c=============================================================
c  main timestepping calculation
c=============================================================
c
c
      do 960 k=1,kmc
      t(k,ic,jc,m,np) = t(k,ic,jc,m,nm)
     &                + c2dtts*boxar*dzr(k)*(
c  u-advection
     &   fuw(k)*(t(k,im,jc,m,nc) + t(k,ic,jc,m,nc))*p5
     & - fue(k)*(t(k,ip,jc,m,nc) + t(k,ic,jc,m,nc))*p5
c  v-advection
     & + fvs(k)*(t(k,ic,jm,m,nc) + t(k,ic,jc,m,nc))*p5
     & - fvn(k)*(t(k,ic,jp,m,nc) + t(k,ic,jc,m,nc))*p5
c  w-advection
     & + fwb1(k  ) - fwb1(k-1)
c  x diffusion
     & +(bbt(jc)*(t(k,ip,jc,m,nm)-t(k,ic,jc,m,nm))*maskpo(k)
     & + bbt(jc)*(t(k,im,jc,m,nm)-t(k,ic,jc,m,nm))*maskmo(k)
c  y diffusion
     & + cct(jc)*(t(k,ic,jp,m,nm)-t(k,ic,jc,m,nm))*maskop(k)
     & + ddt(jc)*(t(k,ic,jm,m,nm)-t(k,ic,jc,m,nm))*maskom(k))
     & * dz(k)
c  z diffusion  
     & + vtf(k  ) - vtf(k-1)
     &   )
 960  continue
 900  continue
c
c=============================================================
c   convectively adjust water column if gravitaionally unstable.
c=============================================================
c
      if(kmc.gt.1.and.ncon.gt.0)then
        do 1070 n=1,ncon
        do 1070 ks=1,2
        call state(t(1,ic,jc,1,np),t(1,ic,jc,2,np),temp(1),
     &              kmc, ks)
        do 1070 m=1,nt
        do 1070 k=ks,kmc-1,2
        if(temp(k).gt.temp(k+1))then
          t(k,ic,jc,m,np) = (dz(k)*t(k,ic,jc,m,np)
     &                    +  dz(k+1)*t(k+1,ic,jc,m,np))*dzw2r(k)
          t(k+1,ic,jc,m,np) = t(k,ic,jc,m,np)
        endif
 1070   continue
      endif
c
      return
      end
 

      subroutine state (t, s, rho, kk, ind)
c
c=============================================================
c
c     state computes one column of normalized densities from the
c     surface to level kk by using a 3rd order polynomial fit to
c     the equation of state.
c
c     note.. for precision purposes, there is a depth dependent
c     constant subtracted from the density returned by this
c     routine.
c     so.. this routine should be used only for calculating  
c     horizontal gradients of density or vertical stability
c     with adjacent levels referred to the same reference level.
c
c     inputs:
c
c     t = the input column of temperatures (degrees C)
c     s = the input column of salinities (units: (ppt-35)/1000)
c     ind = 0 for normal density calculations
c         = 1 for comparing levels 1 to 2, 3 to 4, etc.
c         = 2 for comparing levels 2 to 3, 4 to 5, etc.
c           (if ind = 1 or 2, the coefficients for the lower
c            of the 2 levels are used)
c     kk = bottom level to be calculated
c
c     note: model salinities are equal to (S-35.0)/1000.0 where
c     S is the standard oceanographic salinity.  This convention,
c     which dates from the early Bryan and Cox model, increases 
c     the precision by about 5 binary bits.  This can be 
c     significant in density gradient calculations on 32-bit
c     computers.  The convention is also the one used by the
c     mom/gfdl program used to generate the include file  
c     'dncoef.h'.
c
c     output:
c
c     rho = density with depth level dependent offset 
c           (units: g/cc minus constant)
c
c=============================================================
c
c
#include "param.h"
      dimension t(km),  s(km),  rho(km)
      dimension to(km), so(km), c(km,9),
     &          to1(km),so1(km),co1(km,9),
     &          to2(km),so2(km),co2(km,9)
      save in,to,so,c,to1,so1,co1,to2,so2,co2
#include "dncoef.h"
      data in/1/
c
c  on first entry set up coefficients for stability calculations
c
      if(in.eq.1)then
        in=-1
        do 100 k=1,km-1,2
        to1(k)  = to(k+1)
        to1(k+1)= to(k+1)
        so1(k)  = so(k+1)
        so1(k+1)= so(k+1)
        do 100 l=1,9
        co1(k,l)  = c(k+1,l)
        co1(k+1,l)= c(k+1,l)
 100    continue
        do 200 k=2,km-1,2
        to2(k)  = to(k+1)
        to2(k+1)= to(k+1)
        so2(k)  = so(k+1)
        so2(k+1)= so(k+1)
        do 200 l=1,9
        co2(k,l)  = c(k+1,l)
        co2(k+1,l)= c(k+1,l)
 200    continue
      endif
c
      do 210 k=1,km
        rho(k)=1e30
 210  continue
c
c  test for errors
c
      if(ind.lt.0.or.ind.gt.2.or.kk.lt.0.or.kk.gt.km)then
        write (stderr,99)ind,kk
        stop '1 state'
c
c  normal entry
c
      elseif(ind.eq.0)then
        do 400 k=1,kk
        tq=t(k)-to(k)
        sq=s(k)-so(k)
        rho(k) = (c(k,1) + (c(k,4) + c(k,7)*sq)*sq +
     &           (c(k,3) + c(k,8)*sq + c(k,6)*tq)*tq)*tq +
     &           (c(k,2) + (c(k,5) + c(k,9)*sq)*sq)*sq

 400    continue
c
c  compare 1 with 2, 3 with 4, etc.  density w.r.t bottom level
c
      elseif(ind.eq.1.and.kk.gt.1)then
        do 500 k1=0,kk-2,2
        do 500 k2=1,2
        k=k1+k2
        tq=t(k)-to1(k)
        sq=s(k)-so1(k)
        rho(k) = (co1(k,1) + (co1(k,4)    + co1(k,7)*sq)*sq
     &         + (co1(k,3) +  co1(k,8)*sq + co1(k,6)*tq)*tq)*tq
     &         + (co1(k,2) + (co1(k,5)    + co1(k,9)*sq)*sq)*sq
 500    continue
c
c  compare 2 with 3, 4 with 5 etc
c
      elseif(ind.eq.2.and.kk.gt.2)then
        do 600 k1=1,kk-2,2
        do 600 k2=1,2
        k=k1+k2
        tq=t(k)-to2(k)
        sq=s(k)-so2(k)
        rho(k) = (co2(k,1) + (co2(k,4)    + co2(k,7)*sq)*sq
     &         + (co2(k,3) +  co2(k,8)*sq + co2(k,6)*tq)*tq)*tq
     &         + (co2(k,2) + (co2(k,5)    + co2(k,9)*sq)*sq)*sq
 600    continue
      endif
c
      return
   99 format(/,' Routine state. Variables "ind" or "kk" out',
     &       ' of range',/,' ind =',i10,' kk =',i10)
      end



      subroutine frees(ic,jc)
c
c====================================================================
c
c  Based on the routine in "Killworth, P.D., Stainforth, D.,
c  Webb, D.J. and Paterson, S.M. (1989).  A free surface
c  Bryan-Cox-Semtner model.  Institute of Oceanographic
c  Sciences, Report No. 270.  Wormley, Godlaming, U.K.. 184pp."
c
c=====================================================================
c
c---------------------------------------------------------------------
c  define global data
c---------------------------------------------------------------------
c
#include "param.h"
c
#include "chmix.h"
#include "ctmngr.h"
#include "frees.h"
#include "grdvar.h"
#include "levind.h"
#include "scalar.h"
#include "switch.h"
#include "timelv.h"
c
c
c---------------------------------------------------------------------
c  set indices
c---------------------------------------------------------------------
c
      ip=ic+1
      im=ic-1
      jp=jc+1
      jm=jc-1
c
c---------------------------------------------------------------------
c  timestep height field.
c---------------------------------------------------------------------
c
      if(kmt(ic,jc).ne.0)then
        boxar = dxr*dyr*cstr(jc)
        dhdt  = (dy*(u0(im,jm,nc0)*h(im,jm)+u0(im,jc,nc0)*h(im,jc)
     &              -u0(ic,jm,nc0)*h(ic,jm)-u0(ic,jc,nc0)*h(ic,jc))
     & +dx*(csu(jm)*(v0(im,jm,nc0)*h(im,jm)+v0(ic,jm,nc0)*h(ic,jm))
     &     -csu(jc)*(v0(im,jc,nc0)*h(im,jc)+v0(ic,jc,nc0)*h(ic,jc))
     &     ))*0.5*boxar
c
        h0(ic,jc,np0) = h0(ic,jc,nm0) + dtbt*dhdt
      endif
c
c---------------------------------------------------------------------
c  timestep velocity fields using semi-implicit scheme
c  for the coriolis term.
c---------------------------------------------------------------------
c
      if(kmu(ic,jc).ne.0)then
        fac1=fcor(jc)*dtbt*p5
        fac2=c1/(c1+fac1*fac1)
c
c---------------------------------------------------------------------
c    compute pressure gradients
c---------------------------------------------------------------------
c
        fxa = grav*dx2r*csur(jc)
        fxb = grav*dy2r
        temp1 = h0(ip,jp,nc0) - h0(ic,jc,nc0)
        temp2 = h0(ic,jp,nc0) - h0(ip,jc,nc0)
        dpdx  = (temp1-temp2)*fxa
        dpdy  = (temp1+temp2)*fxb
c
c---------------------------------------------------------------------
c  calculate partial solutions
c---------------------------------------------------------------------
c
        ustar = u0(ic,jc,nm0) + dtbt*(zu(ic,jc) - dpdx)
     &                        + fac1*v0(ic,jc,nm0)
        vstar = v0(ic,jc,nm0) + dtbt*(zv(ic,jc) - dpdy)
     &                        - fac1*u0(ic,jc,nm0)
c
c---------------------------------------------------------------------
c  final step of semi-implicit scheme
c---------------------------------------------------------------------
c
        u0(ic,jc,np0)= (ustar + fac1*vstar)*fac2
        v0(ic,jc,np0)= (vstar - fac1*ustar)*fac2
      endif
c
      return
      end


#ifdef presetp
      subroutine setp (ic, jc)
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
#include "param.h"
#include "grdvar.h"
#include "scalar.h"
#include "slabs.h"
#include "timelv.h"
#include "levind.h"
c
      dimension temp(km)
c
      kmc = kmt(ic,jc)
      if(kmc.gt.0)then        
        call state (t(1,ic,jc,1,nc), t(1,ic,jc,2,nc),  
     &              temp, kmt(ic,jc) ,0)
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
      return
      end
#endif



      subroutine setvbc (ic, jc)
c
c==================================================================
c
c     compute vertical boundary conditions for clinic & tracer
c
c==================================================================
c
#include "param.h"
#include "cvbc.h"
#include "grdvar.h"
#include "scalar.h"
#include "slabs.h"
#include "timelv.h"
#include "levind.h"
c
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


      subroutine bcest (stf, smf, i, j)
c
c==================================================================
c     this routine linearly interpolates global, zonal mean values 
c     of ocean model surface boundary conditions (SST, salinity, 
c     WSX, WSY) having 4.5 deg spacing, to the requested latitdude.
c==================================================================
#include "param.h"
      parameter (nolat = 40, nolatp = nolat+1)
#include "coord.h"
      dimension sstobs(nolat), salobs(nolat)
      dimension wsxobs(nolatp), wsyobs(nolatp)
      dimension olatt(nolat), olatv(nolatp)
      dimension stf(nt), smf(2)
c
c     units of boundary condition values are:
c       degrees C, parts per thousand, dynes per square centimeter
c     olatt  = latitude points for observed data
c     olatv  = latitude points for observed data
c     dolat  = latitude spacing for observed data
c
c     "observed" temperature and salinity data are based on global,
c     annual mean zonally averaged values from the Levitus Atlas
c     (1982).
c     "observed" windstress data are based on global, annual mean,
c     zonally averaged values from Hellerman and Rosenstein (1981).
c     some smoothing was done.
c
c     units:  temperature - degrees Centigrade
c             salinity    - standard oceanographic salinity units
c             wind stress - dynes per sq cm
c     The salinity is transformed into model units in routine
c     'setvbc'.
c
c     references:
c       Hellerman, S, and M. Rosenstein, normal monthly wind stress
c         over the world ocean with error estimates, Journal of
c         Physical Oceanography, 13, 1093-1104,1983.
c       Levitus, S., Climatological atlas of the world ocean, NOAA
c         Prof. Paper 13, US Gov't printing Office, Washington, DC, 
c         1982.
c
      save in,dolat,rdolat,olatt,olatv,sstobs,salobs,wsxobs,wsyobs
c
      data sstobs / -1.75, -1.75, -1.50, -1.50, -1.28,
     &              -0.55,  0.90,  2.92,  5.45,  8.62,
     &              12.27, 15.49, 18.30, 20.67, 22.64,
     &              24.14, 25.27, 26.37, 26.52, 26.16,
     &              26.85, 27.27, 26.82, 26.42, 25.53,
     &              24.03, 22.07, 19.73, 17.02, 12.77,
     &               8.93,  7.25,  6.22,  4.67,  4.57,
     &               3.03, -0.01, -1.05, -1.75, -1.75/
c
      data salobs / 34.30, 34.30, 34.30, 34.13, 33.98,
     &              33.97, 33.97, 33.98, 34.03, 34.24,
     &              34.61, 35.02, 35.37, 35.61, 35.72,
     &              35.68, 35.51, 35.22, 35.05, 35.12,
     &              34.80, 34.56, 34.71, 34.90, 35.27,
     &              35.67, 35.56, 35.49, 35.23, 34.28,
     &              33.57, 33.57, 33.60, 33.80, 34.04,
     &              34.05, 32.65, 32.30, 32.10, 32.00/
c
      data wsxobs /  0.00,
     &               0.00,  0.00, -0.02,  0.15,  0.31,
     &               0.50,  0.82,  1.08,  1.23,  1.16,
     &               0.84,  0.41,  0.02, -0.35, -0.55,
     &              -0.67, -0.64, -0.46, -0.29, -0.19,
     &              -0.16, -0.33, -0.52, -0.59, -0.55,
     &              -0.32,  0.09,  0.42,  0.56,  0.76,
     &               0.81,  0.65,  0.29,  0.06, -0.10,
     &              -0.05, -0.03,  0.05,  0.10,  0.01/
c
      data wsyobs /  .000,
     &               .000,  .009,  .032,  .005, -.023,
     &              -.075, -.155, -.202, -.230, -.179,
     &              -.049,  .093,  .214,  .294,  .344,
     &               .383,  .364,  .269,  .189,  .178,
     &               .125, -.122, -.213, -.251, -.259,
     &              -.202, -.189, -.179, -.183, -.009,
     &               .023,  .053, -.048, -.185, -.225,
     &              -.097, -.050, -.023, -.006,  .000/
c
      data in/1/
c
c----------------------------------------------------------------
c     set latitudes of sst and salinity observations
c     and set latitudes of windstress observations
c----------------------------------------------------------------
c
      if(in.eq.1) then
        in = -1
        dolat = 180.0/nolat
        rdolat = c1/dolat
        do 100 n=1,nolat
          olatt(n) = -90.0 + (n-p5)*dolat
          olatv(n) = -90.0 + (n-1.0)*dolat
100     continue
        olatv(nolatp) = -90.0 + (nolat)*dolat
      endif
c
c----------------------------------------------------------------
c   use linear interpolation to produce the estimated surface 
c   boundary condition values for temperature and salinity at t,s 
c   row j
c----------------------------------------------------------------
c
c
      ylatt = yt(j)
      if (ylatt .le. olatt(1)) then
        nn = 1
        d = c0
      elseif (ylatt .ge. olatt(nolat)) then
        nn = nolat-1
        d = dolat
      else
        do 200 jj=2,nolat
          if (ylatt .le. olatt(jj)) then
            nn = jj - 1
            d  = ylatt - olatt(nn)
            goto 201
          endif
200     continue
      endif
201   continue
      stf(1) = (sstobs(nn)*(dolat - d) + sstobs(nn+1)*d)*rdolat
      stf(2) = (salobs(nn)*(dolat - d) + salobs(nn+1)*d)*rdolat
c
c----------------------------------------------------------------
c   use linear interpolation to produce the estimated surface 
c   boundary condition values for wind stress components at 
c   u,v row j
c----------------------------------------------------------------
c
      ylatv = yu(j)
c
      if (ylatv .le. olatv(1)) then
        nn = 1
        d = c0
      elseif (ylatv .ge. olatv(nolatp)) then
        nn = nolatp - 1
        d = dolat
      else
        do 300 jj=2,nolatp
          if (ylatv .le. olatv(jj)) then
            nn = jj - 1
            d  = ylatv - olatv(nn)
            goto 301
          endif
300     continue
      endif
301   continue
      smf(1) = (wsxobs(nn)*(dolat - d) + wsxobs(nn+1)*d) *rdolat
      smf(2) = (wsyobs(nn)*(dolat - d) + wsyobs(nn+1)*d) *rdolat
c
      return
      end


      subroutine addv( ic, jc)
c
c==================================================================
c
c  add barotropic to baroclinic velocities
c
c==================================================================
c
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


      subroutine diag( ic, jc)
c
c=============================================================
c  collect timestep statistics and (if diagts is set) statistics
c  for diagnostic timestep.
c  note 1. this routine sums arrays, it may not vectorise well
c       2. diagnostic timestep option not implemented yet
c=============================================================
c
c-----------------------------------------------------------------
c         accumulate kinetic energy 
c-----------------------------------------------------------------
c
c
#include "param.h"
c
#include "scalar.h"
#include "switch.h"
#include "grdvar.h"
#include "levind.h"
#include "timelv.h"
#include "slabs.h"
c
#include "cdiag.h"
c
      ku = kmu(ic,jc)
      kt = kmt(ic,jc)
      if(ku.ne.0)then
        boxa = dx*csu(jc)*dy
        do 100 k=1,ku 
        boxv  = boxa*dz(k)
        ektot = ektot + (u(k,ic,jc,nc)**2 + v(k,ic,jc,nc)**2)
     &                  *p5*boxv
 100    continue
      endif
c
c-----------------------------------------------------------------
c         accumulate tracer change statistics 
c-----------------------------------------------------------------
c
      if(kt.ne.0)then
        boxa = cst(jc)*dy*dx/c2dtts
        do 500 m=1,nt
        do 500 k=1,kt
        boxv    = boxa*dz(k)
        tddt(m)= tddt(m)
     &            + (t(k,ic,jc,m,np)-t(k,ic,jc,m,nm))*boxv
        dtabs(m)= dtabs(m)
     &            + abs(t(k,ic,jc,m,np)-t(k,ic,jc,m,nm))*boxv
        tvar(m) = tvar(m)+(t(k,ic,jc,m,np)**2
     &                       -t(k,ic,jc,m,nm)**2)*boxv
500     continue
      endif
      return
      end

      subroutine grids
c==================================================================
c
c     calculate all grid related items
c
c==================================================================
c
#include "param.h"
#include "scalar.h"
#include "coord.h"
#include "grdvar.h"
c
c------------------------------------------------------------------
c     set latitudes & longitudes of "t" & "u,v" grid points (in degrees)
c------------------------------------------------------------------
c
      xt(1) = stlon + p5*dxdeg
      xu(1) = stlon + dxdeg
      do 100 i=2,imt
        xt(i) = xt(i-1) +  dxdeg
        xu(i) = xu(i-1) +  dxdeg
100   continue
c
      yt(1) = stlat + p5*dydeg
      yu(1) = stlat + dydeg
      do 200 j=2,jmt
        yt(j) = yt(j-1) +  dydeg
        yu(j) = yu(j-1) +  dydeg
200   continue
c
c     convert to centimeters
c
      dy = dydeg*radius/radian
      dx = dxdeg*radius/radian
c
c------------------------------------------------------------------
c     compute auxiliary arrays
c------------------------------------------------------------------
c
      do 500 k=1,km
        dzr(k)  = c1/dz(k)
        c2dz(k) = c2*dz(k)
        dz2r(k) = c1/c2dz(k)
500   continue
c
      dzw(0) = p5*dz(1)
      zw(1)  = dz(1)
      do 600 k=2,km
        dzw(k-1) = p5*(dz(k-1)+dz(k))
        zw(k)    = zw(k-1)+dz(k)
600   continue
c
      dzw(km)   = p5*dz(km)
      dzwr(km)  = c1/dzw(km)
      dzw2r(km) = p5/dzw(km)
      zt(1)     = dzw(0)
c
      do 700 k=1,km
        dzwr(k-1)  = c1/dzw(k-1)
        dzw2r(k-1) = p5/dzw(k-1)
        zt(k+1)  = zt(k) + dzw(k)
700   continue
c
      dxr = c1/dx
      dx2r = p5/dx
      dx4r = p25/dx
c
      dyr = c1/dy
      dy2r = p5/dy
      dy4r = p25/dy
c
      phi(1)   = yu(1)/radian
      phit(1)  = phi(1) - p5*dydeg/radian
      sumdy    = phi(1)
c
      do 900 j=1,jmt
        if (j .ne. jmt) sumdy = sumdy + dydeg/radian
        if (j .ne. jmt) phi(j+1) = sumdy
        if (j .ne. 1) phit(j) = p5*(phi(j-1)+phi(j))
        cst(j)  = cos(phit(j))
        csu(j)  = cos(phi (j))
        sine(j) = sin(phi(j))
        cstr(j) = c1/cst(j)
        csur(j) = c1/csu(j)
        tng(j)  = sine(j)/csu(j)
        tanra(j)= tng(j)/radius
900   continue
c
c------------------------------------------------------------------
c     print grid related arrays
c------------------------------------------------------------------
c
c
      return
c
9001  format (1x,10g13.5)
9101  format (/,' "t" grid box thickness (cm): "dz(k) k=1..km"')
9102  format (/,' "w" grid box thickness (cm): ',
     &        '"dzw(k) k=0..km"')
9103  format (/,' depth to "t" box bottom (cm): "zw(k) k=1..km"')
9104  format (/,' depth to "t" & "u,v" grid points (cm): "zt(k) ',
     &        'k=1..kmp1"')
9105  format (/,' latitude of t,s points (deg): "yt(j) j=1,jmt"')
9106  format (/,' latitude of u,v points (deg): "yu(j) j=1,jmt"')
9107  format (/,' longitude of t,s points (deg): "xt(i) i=1,imt"')
9108  format (/,' longitude of u,v points (deg): "xu(i) i=1,imt"')
      end


      subroutine ocn1st
c
c==================================================================
c
c  initialize quantites when model first started
c
c==================================================================
c
#include "param.h"
#include "coord.h"
#include "ctmngr.h"
#include "cvbc.h"
#include "frees.h"
#include "grdvar.h"
#include "iounit.h"
#include "levind.h"
#include "scalar.h"
#include "slabs.h"
#include "switch.h"
c
c------------------------------------------------------------------
c     generate topography over "t" points
c------------------------------------------------------------------
c
      call topog
c
c------------------------------------------------------------------
c     initialize t,s,u,v.
c          velocities set to zero
c          temperature an salinity set to 0.0 and 35.0 at sea
c                                      to -2.0 and 45.0 over land.
c------------------------------------------------------------------
c
      sconst=0.0349
      do 200 n=1,3
      do 200 j=1,jmt
      do 200 i=1,imt
      do 200 k=1,km
        u(k,i,j,n)=c0
        v(k,i,j,n)=c0
        if(k.gt.kmt(i,j))then
          t(k,i,j,1,n)=-2.0
          t(k,i,j,2,n)=0.01
        else
          t(k,i,j,1,n)= theta0 (yt(j), zt(k))
          t(k,i,j,2,n)= sconst-0.035
      endif
200   continue
c
c------------------------------------------------------------------
c     initialize free surface model
c------------------------------------------------------------------
c
      do 300 n=1,3
      do 300 j=1,jmt
      do 300 i=1,imt
        h0(i,j,n) = c0
        u0(i,j,n) = c0
        v0(i,j,n) = c0
300   continue
c
c------------------------------------------------------------------
c     initialize controls
c------------------------------------------------------------------
c
      itt   = 0
      totsec = c0
      totday = c0
      years  = c0
c
      return
      end


      subroutine topog
c
c------------------------------------------------------------------
c     construct an "idealized" world ... piece by piece
c     note: this topography will map into arbitrary model 
c     resolution
c------------------------------------------------------------------
c
#include "param.h"
#include "coord.h"
#include "grdvar.h"
#include "levind.h"
c
c     define the model topography on the "t" grid using "kmt"
c     set maximum level everywhere except on boundary
c
      do 100 j=2,jmtm1
        do 90 i=2,imtm1
          kmt(i,j) = km
90      continue
100   continue
c
      do 200 j=1,jmt
        kmt(imt,j) = 0
        kmt(1,j) = 0
200   continue
      do 300 i=1,imt
        kmt(i,jmt) = 0
        kmt(i,1) = 0
300   continue
c
c     antarctica
c
      call setkmt (-90.0, 0.0, c360, -80.0, 0.0, c360, 0)
      call setkmt (-80.0, c360-25.0, c360, -70.0, c360, c360, 0)
      call setkmt (-80.0, 0.0, c360, -70.0, 0.0, 170.0, 0)
      call setkmt (-80.0, c360-135.0, c360-60.0, -68.0, c360-75.0,
     &              c360-60.0, 0)
      call setkmt (-70.0, 0.0, 155.0, -67.0, 50.0, 145.0, 0)
c
c     australia
c
      call setkmt (-35.0, 116.0, 120.0, -31.0, 114.0, 130.0, 0)
      call setkmt (-38.0, 140.0, 151.0, -31.0, 130.0, 151.0, 0)
      call setkmt (-31.0, 115.0, 153.0, -20.0, 113.0, 149.0, 0)
      call setkmt (-20.0, 113.0, 149.0, -11.0, 131.0, 143.0, 0)
c
c     south america
c
      call setkmt (-50.0, c360-74.0, c360-68.0, -40.0, c360-73.0,
     &              c360-62.0, 0)
      call setkmt (-40.0, c360-73.0, c360-62.0, -20.0, c360-70.0,
     &              c360-40.0, 0)
      call setkmt (-20.0, c360-70.0, c360-40.0, -16.0, c360-81.0,
     &              c360-35.0, 0)
      call setkmt (-16.0, c360-81.0, c360-35.0, 0.0, c360-80.0,
     &             c360-50.0, 0)
      call setkmt (0.0, c360-80.0, c360-50.0, 11.0, c360-75.0,
     &             c360-60.0, 0)
c
c     central america
c
      call setkmt (6.0, c360-78.0, c360-75.0, 20.0, c360-105.0,
     &             c360-97.0, 0)
      call setkmt (20.0, c360-105.0, c360-97.0, 30.0, c360-115.0,
     &             c360-94.0, 0)
c
c     north america
c
      call setkmt (25.0, c360-82.0, c360-80.0, 30.0, c360-85.0,
     &             c360-81.0, 0)
      call setkmt (30.0, c360-115.0, c360-80.0, 40.0, c360-124.0,
     &             c360-74.0, 0)
      call setkmt (40.0, c360-124.0, c360-74.0, 50.0, c360-124.0,
     &             c360-57.0, 0)
      call setkmt (50.0, c360-124.0, c360-57.0, 60.0, c360-140.0,
     &             c360-64.0, 0)
      call setkmt (60.0, c360-165.0, c360-64.0, 65.0, c360-140.0,
     &             c360-64.0, 0)
      call setkmt (65.0, c360-140.0, c360-64.0, 70.0, c360-162.0,
     &             c360-72.0, 0)
      call setkmt (70.0, c360-162.0, c360-140.0, 72.0, c360-157.0,
     &             c360-157.0, 0)
      call setkmt (70.0, c360-130.0, c360-70.0, 75.0, c360-120.0,
     &             c360-80.0, 0)
c
c     greenland
c
      call setkmt (60.0, c360-45.0, c360-45.0, 75.0, c360-58.0,
     &             c360-19.0, 0)
c
c     africa
c
      call setkmt (-35.0, 19.0, 28.0, 6.0, 8.0, 50.0, 0)
      call setkmt (6.0, 0.0, 50.0, 18.0, 0.0, 56.0, 0)
      call setkmt (18.0, 0.0, 56.0, 26.0, 0.0, 59.0, 0)
      call setkmt (6.0, c360-10.0, c360, 18.0, c360-18.0, c360, 0)
      call setkmt (18.0, c360-18.0, c360, 26.0, c360-15.0, c360, 0)
c
c     northern africa &  europe & asia
c
      call setkmt (26.0, c360-15.0, c360, 40.0, c360-7.0, c360, 0)
      call setkmt (40.0, c360-7.0, c360, 50.0, c360, c360, 0)
c
      call setkmt (26.0, 0.0, 126.0, 40.0, 0.0, 122.0, 0)
      call setkmt (40.0, 0.0, 130.0, 50.0, 0.0, 140.0, 0)
      call setkmt (50.0, 0.0, 140.0, 60.0, 8.0, 140.0, 0)
      call setkmt (60.0, 8.0, 163.0, 65.0, 13.0, 180.0, 0)
      call setkmt (65.0, 13.0, 188.0, 70.0, 20.0, 180.0, 0)
      call setkmt (70.0, 70.0, 180.0, 75.0, 90.0, 100.0, 0)
c
c     add an "idealized" ridge in the atlantic & pacific
c
      level = indp (3500.e2, zw, km)
c
      call setkmt (-20.0, c360-20.0, c360-10.0, 30.0, c360-45.0,
     &              c360-35.0, level)
      call setkmt (30.0, c360-45.0, c360-35.0, 60.0, c360-20.0,
     &              c360-30.0, level)
c
      call setkmt (-60.0,c360-100.0, c360-130.0, 40.0, c360-160.0,
     &              180.0, level)
      level = indp (2000.e2, zw, km)
      call setkmt (-50.0, c360-120.0, c360-120.0, 30.0, 190.0, 190.0,
     &             level)
c
c     set cyclic boundary conditions
c
      do 400 j=1,jmt
        kmt(1,j)   = kmt(imtm1,j)
        kmt(imt,j) = kmt(2,j)
400   continue
      return
      end


      function indp(c, ca, ia)
c
c
c==================================================================
c
c     indp = index of nearest data point
c
c     inputs:
c
c     c      = arbitrary data point...same units as ca
c     ca     = array of data points
c              (must be monotonically increasing)
c     ia     = dimension of ca
c
c     output:
c
c     indp =  index of nearest data point to "c"
c             if "c" is outside the domain of "ca" then i
c             ndp is set to 1 or ia depending on  whether
c             c(1) or c(ia) is closest
c
c        note: if "ca" is dimensioned ca(0:ia) in the calling
c              program, then the returned index should be reduced
c              by one to account for the zero base.
c
c     example:
c
c     let model depths be defined by the following:
c     parameter (km=5)
c     dimension z(km)
c     data z /5.0, 10.0, 50.0, 100.0, 250.0/
c
c     k1 = indp (12.5, z, km)
c     k2 = indp (0.0, z, km)
c
c     k1 would be set to 2, & k2 would be set to 1 so that
c     z(k1) would be the nearest data point to 12.5 and z(k2) would
c     be the nearest data point to 0.0
c
c==================================================================
c
c
      dimension ca(ia)
c
      if (c .lt. ca(1) .or. c .gt. ca(ia)) then
        if (c .lt. ca(1))  indp = 1
        if (c .gt. ca(ia)) indp = ia
        return
      else
        do 10 i=2,ia
          if (c .le. ca(i)) then
            indp = i
            if (ca(i)-c .gt. c-ca(i-1)) indp = i-1
            go to 11
          endif
10      continue
11      continue
      endif
      return
      end


      subroutine setkmt (alat1, slon1, elon1,
     &                   alat2, slon2, elon2, num)
c
c------------------------------------------------------------------
c     set the topography mask "kmt(i,j)" = "num" within the area of 
c     the parallelogram bounded by vertices:
c     (alat1,slon1), (alat1,elon1), (alat2,slon1), & (alat2,elon2)
c------------------------------------------------------------------
c
#include "param.h"
#include "coord.h"
#include "levind.h"
c
c     convert the four vertices into model indices
c     (js,is1), (js,ie1), (je,is2), (je,ie2)
c
      j1 = indp (alat1, yt, jmt)
      j2 = indp (alat2, yt, jmt)
      js = min (j1,j2)
      je = max (j1,j2)
c
      i1  = indp (slon1, xt, imt)
      i2  = indp (elon1, xt, imt)
      is1 = min (i1,i2)
      ie1 = max (i1,i2)
c
      i1  = indp (slon2, xt, imt)
      i2  = indp (elon2, xt, imt)
      is2 = min (i1,i2)
      ie2 = max (i1,i2)
c
      is = is1
      ie = ie1
c
c     fill in the area bounded by (js,is1), (js,ie1), 
c                                 (je,is2), (je,ie2)
c
      if (js .eq. je) then
        rdj = c1
      else
        rdj = c1/(je-js)
      endif
      do 100 j=js,je
        do 90 i=is,ie
          kmt(i,j) = num
90      continue
        is = nint(rdj*((j-js)*is2 + (je-j)*is1))
        ie = nint(rdj*((j-js)*ie2 + (je-j)*ie1))
100   continue
      return
      end


      function theta0 (ydeg, depth)
c
c==================================================================
c     this subroutine returns estimates of global mean potential
c     temperature for model initialization as a function of depth.
c     it is used to produce a reference thermal stratification for 
c     the upper 2000m of the MOM's test case.  below 2000m, the
c     potential temperature returned is 2.0 degrees C.  surface
c     values are set slightly above 18.4 degrees C at the reference
c     latitude "reflat".
c     the estimates are produced from a 7th order ploynomial fit to
c     the annual mean world ocean potential temperature observations
c     of Levitus (1982).
c
c     input [units]:
c       a latitdue (ydeg): [degrees]
c       a zt value (depth): [centimeters]
c     output [units]:
c       potential temperature estimate (est): [degrees centigrade]
c
c     variables:
c       coeft     = coefficients for the polynomial fit of potential
c                   temperature vs. depth
c       reflat    = reference latitude at which observed surface
c                   temperatures approximately equal coeft(1)
c       factor    = the ratio of the cosine of the latitude 
c                   requested ("ydeg") to the reference latitude 
c                   ("reflat")used to scale the upper 2000 meters
c                   of the vertical temperature profile
c       tmin,tmax = the minumum and maximum potential temperatures
c                   allowed at the time of model initialization
c
c     reference:
c       Levitus, S., Climatological atlas of the world ocean, NOAA
c     Prof. Paper 13, US Gov't printing Office, Washington, DC,
c     1982.
c
      parameter (ndeg=7)
      dimension coeft(ndeg+1)
      save coeft, tmin, tmax, reflat, c0,c1,c2,pi,refcos
      data coeft / 0.184231944E+02,-0.430306621E-01, 0.607121504E-04
     &           ,-0.523806281E-07, 0.272989082E-10,-0.833224666E-14
     &           , 0.136974583E-17,-0.935923382E-22/
      data tmin, tmax, reflat /2.0, 25.0, 34.0/
      data in/1/
c
c==================================================================
c
      if(in.eq.1)then
      in=-1
      c0 = 0.0
      c1 = 1.0
      c2 = 2.0
      pi = atan(1.0) * 4.0
      refcos = abs(cos(pi*reflat/180.))
      endif
c
      coslat = abs(cos(pi*ydeg/180.))
      factor = coslat/refcos
      z = depth * 0.01
c
      if (z .gt. 2000.) then
        est = c2
      else
        est = c0
        zz  = c1
        do 100 nn=1,ndeg+1
          est = est + coeft(nn)*zz
          zz = zz*z
100     continue
        est = est * factor
      endif
c
      if (est .gt. tmax) est = tmax
      if (est .lt. tmin) est = tmin
c
      theta0 = est
      return
      end


      subroutine tmngr (dt)
c
c
c==================================================================
c
c                      t-i-m-e      m-a-n-a-g-e-r
c
c     inputs:
c
c     itt    = current model time step (from ctmngr.h).
c     dtts   = length of time step in seconds.
c     dgnstc = number of days between printing diagnostics
c     nmix   = number of timesteps between mixint timesteps
c     tsi    = number of days between printing timestep information
c     days   = day at which to end run.
c
c     outputs:
c
c     itt    = itt (on input)  + 1
c     totday = total number of days elapsed 
c              (0.5 = noon on first day).
c     years is number of years elapsed (365.5 days = 1 year)
c     diagts = true if time to print diagnostics
c     mixts  = true if a mixing timestep
c     prntsi = true if time to print timestep information
c     last   = true if last timestep of run
c
c==================================================================
c
#include "pconst.h"
#include "ctmngr.h"
#include "switch.h"
#include "timelv.h"
      logical tmnset
c
c------------------------------------------------------------------
c     increment time step counters
c     and create a time stamp
c------------------------------------------------------------------
c
      itt = itt + 1
      totsec = totsec + dt
      totday = totsec*secday
      years  = totday/365.25
      write(stamp,11)totday,years
  11  format('totday =',f8.2,',years =',f8.3)
c
c------------------------------------------------------------------
c
c     set time dependent logical switches.
c     1. set last if it is the last timestep
c     2. set diagts if time to do diagnostics 
c     3. set mixts if time mixing to be done
c     4. set prntsi if time to print timestep info
c     5. set snapts if time to archive model
c
c------------------------------------------------------------------
c
      last   = (totday .ge. days) 
      diagts = tmnset (dgnstc, totday, dt)
      mixts  = (mod(itt,nmix) .eq. 1) 
      prntsi = tmnset (tsi, totday, dt)
      snapts = tmnset (snaps, totday, dt) .and. .not. last
c
      return
      end


      function tmnset (days, totday, dt)
c
c==================================================================
c     tmnset = true if current time is within 1/2 time step of
c              desired interval
c     days   = desired time interval 
c     totday = accumulated time in days from start of run
c     dt     = time step in seconds
c==================================================================
c
#include "pconst.h"
      logical tmnset
      close  = p5*dt
      i      = nint (totday/days)
      remain = (totday - i*days)*c60*c1440
      if ((remain .ge. c0 .and. remain .lt. close) .or.
     &    (remain .lt. c0 .and. -remain .le. close)) then
        tmnset = .true.
      else
        tmnset = .false.
      endif
      return
      end


      subroutine arch
c
c==================================================================
c     write archive file
c     NOTE - this version writes a standard restart file only
c==================================================================
c
      call wrrest
      return
      end


      subroutine pdiag
c
c==================================================================
c     write diagnostics
c==================================================================
c
#include "param.h"
#include "coord.h"
#include "ctmngr.h"
#include "frees.h"
#include "timelv.h"
#include "slabs.h"
      dimension ptemp(imt,km)
      character*12 trname(nt)
c
c==================================================================
c based on MOM common "cprnts.h" - only used by this routine
c
      parameter (nlatpr=10)
c
      dimension prlat(nlatpr), prslon(nlatpr), prelon(nlatpr)
     &,         prsdpt(nlatpr), predpt(nlatpr)
     &,         slatxy, elatxy, slonxy, elonxy
c
c     default starting & ending coordinates for matrix prints during
c     diagnostic time steps (lats & lons in degrees, depths in cm)
c     see "cprnts.h"
c
      data prlat /-60.0, 0.0, 27.0, 55.0, 6*100.0/
      data prslon /180.0, 180.0, 180.0, 180.0, 6*0.0/
      data prelon /250.0, 250.0, 250.0, 250.0, 6*0.0/
      data prsdpt /nlatpr*0.0/
      data predpt /nlatpr*6000.e2/
      data slatxy, elatxy, slonxy, elonxy /-90.0, 90.0, 3.0, 357.0/
c==================================================================
      data trname/'temperature ','salinity    '/
c
c------------------------------------------------------------------
c       print the solution at specified latitudes
c------------------------------------------------------------------
c
      do 700 jlat=1,nlatpr
      jj = indp (prlat(jlat), yt, jmt)
      if (prlat(jlat) .le. yt(jmt)) then
        if (jlat .eq. 1) write(stdout,9000)
        is = indp (prslon(jlat), xt, imt)
        ie = indp (prelon(jlat), xt, imt)
        ks = indp (prsdpt(jlat), zt, km)
        ke = indp (predpt(jlat), zt, km)
        fx = 1.0e-2
        do 680 m=1,nt
        scl = c1
        if (m .eq. 2) scl=1.e-3
        do 710 k=1,km
        do 710 i=1,imt
 710    ptemp(i,k) = t(k,i,jj,m,nc)
        write (stdout,9100) trname(m), itt
     &,    jj, yt(jj), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
        call matrix (ptemp(1,1), imt, is, ie, ks, ke, scl)
680     continue
c
        do 720 k=1,km
        do 720 i=1,imt
 720    ptemp(i,k) = u(k,i,jj,nc)
        scl = c1
        write (stdout,9100) 'u velocity', itt
     &,   jj, yt(jj), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
        call matrix (ptemp(1,1), imt, is, ie, ks, ke, scl)
c
        do 730 k=1,km
        do 730 i=1,imt
 730    ptemp(i,k) = v(k,i,jj,nc)
        scl = c1
        write (stdout,9100) 'v velocity', itt
     &,   jj, yt(jj), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
        call matrix (ptemp(1,1), imt, is, ie, ks, ke, scl)
      endif
700   continue
c
c-------------------------------------------------------------------
c       print the free surface model solution
c-------------------------------------------------------------------
c
      scl = c100
      is = indp (slonxy, xt, imt)
      ie = indp (elonxy, xt, imt)
      js = indp (slatxy, yt, jmt)
      je = indp (elatxy, yt, jmt)
      write (stdout,8000) ' surface height (m)'
     &, itt, xt(is), xt(ie), yt(js), yt(je), scl
      call matrix (h0(1,1,np0), imt, is, ie, -js, -je, scl)
c
      scl = c1
      write (stdout,8000) ' u velocity (cm/s)'
     &, itt, xt(is), xt(ie), yt(js), yt(je), scl
      call matrix (u0(1,1,np0), imt, is, ie, -js, -je, scl)
c
      write (stdout,8000) ' v velocity (cm/s)'
     &, itt, xt(is), xt(ie), yt(js), yt(je), scl
      call matrix (v0(1,1,np0), imt, is, ie, -js, -je, scl)
c
      return
c
8000  format(1x,a27,1x,'ts=',i7
     &,', lon:',f6.2,' ==> ',f6.2,',   lat:',f6.2,' ==> ',f6.2
     &,', scaling=',1pg10.3)
9000  format(/' Zonal section printouts at specified latitudes:'/)
9100  format(1x,a27,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     &,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     &,', scaling=',1pg10.3)
      end


      subroutine restio
c
c==================================================================
c
c     read and write restart file
c     NOTE - this version only saves and restores the current
c     timestep.  The model must therefore be restarted with
c     a forward or euler backward timestep.
c==================================================================
c
#include "param.h"
#include "ctmngr.h"
#include "frees.h"
#include "iounit.h"
#include "levind.h"
#include "slabs.h"
#include "switch.h"
#include "timelv.h"
      dimension t1(imt,jmt,km,nt), t2(imt,jmt,km,nt), 
     &          t3(imt,jmt,km,nt)
      dimension u1(imt,jmt,km),  u2(imt,jmt,km), u3(imt,jmt,km)
      dimension v1(imt,jmt,km),  v2(imt,jmt,km), v3(imt,jmt,km)
      dimension h01(imt,jmt),    h02(imt,jmt),   h03(imt,jmt)
      dimension u01(imt,jmt),    u02(imt,jmt),   u03(imt,jmt)
      dimension v01(imt,jmt),    v02(imt,jmt),   v03(imt,jmt)
      character*80 file
c
      equivalence (t1(1,1,1,1),t(1,1,1,1,1))
      equivalence (t2(1,1,1,1),t(1,1,1,1,2))
      equivalence (t3(1,1,1,1),t(1,1,1,1,3))
      equivalence (u1(1,1,1),u(1,1,1,1))
      equivalence (u2(1,1,1),u(1,1,1,2))
      equivalence (u3(1,1,1),u(1,1,1,3))
      equivalence (v1(1,1,1),v(1,1,1,1))
      equivalence (v2(1,1,1),v(1,1,1,2))
      equivalence (v3(1,1,1),v(1,1,1,3))
c
      equivalence (h01(1,1),h0(1,1,1))
      equivalence (h02(1,1),h0(1,1,2))
      equivalence (h03(1,1),h0(1,1,3))
      equivalence (u01(1,1),u0(1,1,1))
      equivalence (u02(1,1),u0(1,1,2))
      equivalence (u03(1,1),u0(1,1,3))
      equivalence (v01(1,1),v0(1,1,1))
      equivalence (v02(1,1),v0(1,1,2))
      equivalence (v03(1,1),v0(1,1,3))
c
c------------------------------------------------------------------
c     read restart data set from archive for timesteps itt 
c------------------------------------------------------------------
c
      entry rdrest
c  check np equals 3
      if(np.ne.3)then
        write(stdout,21) np, nc, nm
        stop
      endif
c
c------------------------------------------------------------------
c     initialize t,s,u,v.
c          velocities set to zero
c          temperature and salinity set to land values
c          -2.0 and 45.0 ppt
c------------------------------------------------------------------
c
      do 200 n=1,3
      do 200 j=1,jmt
      do 200 i=1,imt
      do 200 k=1,km
        u(k,i,j,n)=c0
        v(k,i,j,n)=c0
        t(k,i,j,1,n)= -2.0
        t(k,i,j,2,n)= 0.01
200   continue
c
c------------------------------------------------------------------
c     initialize free surface model
c------------------------------------------------------------------
c
      do 300 n=1,3
      do 300 j=1,jmt
      do 300 i=1,imt
        h0(i,j,n) = c0
        u0(i,j,n) = c0
        v0(i,j,n) = c0
300   continue
c
c------------------------------------------------------------------
c  read archived data.  Use file fnrest if defined.
c------------------------------------------------------------------
c
      if(fnrest(1:1).ne.' ')then
      open(iorest, file=fnrest, form='UNFORMATTED', 
     &       access='SEQUENTIAL')
      else
        open(iorest, form='UNFORMATTED', access='SEQUENTIAL')
      endif
      read(iorest) itt, totsec, totday, years, kmt
      read(iorest) h03,  u03,  v03
      read(iorest) t3,   u3,   v3
      close(iorest)
      write (stdout,22) iorest, itt
      if(mod(itt,nmix).ne.0)then
        write(stdout,23) itt, nmix
        stop
      endif
      return
c
  21  format(' Routine restio.',/,
     &       ' Trying to restart with np not equal to 3.',/,
     &       ' np, nc, nm  =',3i5,/,' Programme stopped.')
  22  format(' Routine restio',/,
     &       ' Restart file read from unit ',i2,'.  timestep =',i7)
  23  format(' Routine restio - RESTART ERROR!',/,
     &       ' Trying to restart without forward timestep,',/,
     &       ' i.e. mod(itt,nmix).ne.0.    itt, nmix =',2i5,/,
     &       ' Programme stopped.')
c
c------------------------------------------------------------------
c     write restart data set to archive.
c------------------------------------------------------------------
c
      entry wrrest
c
c------------------------------------------------------------------
c  create archive file name
c------------------------------------------------------------------
c
      write (file,'(i10.10)')nint(totday)
      file = 'd'//file(7:10)//'.restart'
      open(iorest, file=file, form='UNFORMATTED',
     &       access='SEQUENTIAL')
c
      write (stdout,24) itt, file, iorest
      write(iorest) itt, totsec, totday, years, kmt
      if(np0.eq.1) write(iorest)  h01, u01, v01
      if(np0.eq.2) write(iorest)  h02, u02, v02
      if(np0.eq.3) write(iorest)  h03, u03, v03
c
      if(np.eq.1) write(iorest)  t1, u1, v1
      if(np.eq.2) write(iorest)  t2, u2, v2
      if(np.eq.3) write(iorest)  t3, u3, v3
c
      close(iorest, status = 'KEEP')
      return
c
  24  format('  Writing restart file for itt = ',i7,/,
     &       '  filename = ',a,'  unit = ',i2) 
      end


      subroutine matrix (array, irdim, istrt, im, jstrt, jm, scale)
c
c==================================================================
c
c     matrix is a general two-dimensional array printing routine,
c     where:
c     array = the array to be printed
c     irdim = the 1st dimension of array
c     istrt = the 1st element of the 1st dimension to be printed
c     im    = the last element of the 1st dimension to be printed
c     jstrt = the 1st element of the 2nd dimension to be printed
c     jm    = the last element of the 2nd dimension to be printed
c             the 2nd dimension is printed in reverse order if both
c             jstrt & jm are negative
c     scale = a scaling factor by which array is divided before
c             printing.  (if this is zero, no scaling is done.)
c             if scale=0, 10 columns are printed across in e format
c             if scale>0, 20 columns are printed across in f format
c
c==================================================================
c
#include "param.h"
c
      dimension array(irdim,1000)
c
      if (jstrt*jm .lt. 0) then
        write (stderr,999)  jstrt, jm
        stop '=>matrix'
      endif
c
c     allow for inversion of 2nd dimension
c
      if (jm .lt. 0) then
        js   = -jm
        je   = -jstrt
        jinc = -1
      else
        js   = jstrt
        je   = jm
        jinc = 1
      endif
c
      if (scale .eq. c0) then
        do 100 is=istrt,im,10
          ie = min(is + 9,im)
          write (stdout,9001) (i, i=is,ie)
          do 90 l=js,je,jinc
            write (stdout,9002) l, (array(i,l),i=is,ie)
90        continue
          write (stdout,'(//)')
100     continue
c
      else
        scaler = c1/scale
        do 200 is=istrt,im,20
          ie = min(is + 19,im)
          write (stdout,9003) (i, i=is,ie)
          do 190 l=js,je,jinc
            write (stdout,9004) l, (array(i,l)*scaler,i=is,ie)
190       continue
          write (stdout,'(//)')
200     continue
      endif
c
      return
c
999   format (1x,'jstrt=',i5,' jm=',i5,' in matrix')
9001  format(10i13)
9002  format(1x,i2,10(1pe13.5))
9003  format(3x,20i6)
9004  format(1x,i3,1x,20f6.2)
c
      end
