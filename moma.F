      program moma
!
!======================================================================
!
!    File:  moma.F
!    Copyright 1992-1997, 2020  David J. Webb
!
!    This file is part of the Moma ocean model code.
!    Released under licence GPL-3.0-or-later
!
!    Moma is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Moma is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Moma.  If not, see <https://www.gnu.org/licenses/>.
!
!======================================================================
!
!  Program MOMA,     Version 1.17
!
!  Free surface ocean model for use with array processors.
!
!  This code is a free surface version of the Bryan-Cox-Semtner
!  code reorganised for efficient running on array processor
!  computers.  The basic finite difference scheme used in the
!  model is described in Bryan (1969), Semtner (1974) and
!  Cox (1984).  The present code follows most closely
!  the format of the GFDL Modular Ocean Model (Pacanowski et al,
!  1990), and makes use of many of the subroutines and include
!  files of that model.  More detailed information about the
!  present code is given in Webb (1993, 1996).
!
!  The code may be freely adapted and used without charge.  For
!  further information on the model contact the author.
!
!     Copyright 1993, 1994, 1995, 1996
!     D.J. Webb,
!     Southampton Oceanography Centre
!     Empress Dock, Southampton SO14 3ZH, U.K..
!
!     e-mail  David.Webb@soc.soton.ac.uk
!
!  The main differences from the modular ocean model code are:
!
!  1.  Collection of all 'array processor' loops (i.e. loops over the
!      horizontal indices ic and jc) in subroutine step.  When an array
!      processor is used, code should be added to this subroutine to
!      partition the horizontal index (ic, jc) ranges between the
!      different processors so that they each have a similar workload.
!  2.  Revised common block structures.  All variables are stored
!      in core.
!  3.  Removal of all 'slab' optimisation code in the program and
!      the introduction of code designed to optimise inner loop
!      calculations using the vertical index k.
!  4.  To simplify the development and testing of the array processor
!      code the following features of the moma code were removed:
!        (a) Diagnostic calculations in routines clinic and tracer.
!        (b) Moma code options.
!      It should now be straightforward (but possibly time consuming)
!      to add the features back into the present code if they are
!      required:
!  5.  The stream function code is removed and replaced by a free-
!      surface code.  The latter is similar to Killworth, Stainforth
!      Webb and Patterson (1989).  For efficiency, the mean horizontal
!      velocity is now used instead of horizontal transport and
!      the viscous terms are calculated in routine clinic.
!      A leapfrog timestepping scheme is used with averaging over each
!      baroclinic timestep to prevent aliasing problems in the tropics.
!  6.  In the baroclinic momentum equation, a revised horizontal
!      advection scheme option is included.  The old scheme can still be
!      used by specifying the flag 'oldadv'
!  7.  Near topography, in the baroclinic momentum equation, a revised
!      vertical advection scheme option is included.
!  8.  An option to precalculate the baroclinic part of the
!      pressure filed is included using flag 'presetp'.
!  9.  Snapshot and archive functions separated.  Input 'snapd' and 'archd'.
!
!
!  cpp precompiler options:
!
!     'oldadv'    - use the origonal scheme for the horizontal
!                   advection velocity at velocity points.
!     'presetp'   - precalculate the baroclinic pressure field
!                   before calling clinic
!     'hcomments' - include comments from '*.h' files
!     'cyclic'    - invoke cyclic boundary conditions
!     'simple_world' - use routine topog to construct kmt, this generates
!                   a simplified version of the global ocean.
!     'free_eb'   - use the euler backward timestepping scheme without
!                   averaging for the free surface equations.
!                   NOTE: The default is to use a leapfrog scheme
!                   with time averaging.
!     'de_checkbd'- use filter to remove checkerboarding in free
!                   surface height field.
!     'ncsa_hdf'  - produce hdf versions of archive files.
!
!  Compilation:
!
!      Without 'ncsa_hdf' the program can be compiled with a command
!      of the form
!             f77 -o moma moma.F
!      This assumes that the fortran compiler also calls the cpp
!      preprocessor.  With 'ncsa_hdf' the hdf include files and
!      subroutines are also needed.  These require commands like:
!
!             set HDFROOT = /users/occam/hdf/sunos2
!             f77 -Dncsa_hdf -I$HDFROOT/include -o moma moma.F
!                  $HDFROOT/lib/libnetcdf.a $HDFROOT/lib/libdf.a
!
!      where the third line is really a continuation of the f77 command line
!      and HDFROOT is set to point to the hdf root directory.
!
!  References:
!
!  Bryan, K., 1969:  A numerical method for the circulation of the
!  World Ocean.  Journal of Computational Physics, 4, 347-  .
!
!  Semtner, A.J., 1974:  A general circulation model for the
!  World Ocean.  UCLA Department of Meteorology Technical Report
!  No. 8, 99pp.
!
!  Cox, M.D., 1984:  A primitive equation, 3-dimensional model of
!  the ocean.  GFDL Ocean Technical Report No.1, Geophysical
!  Fluid Dynamics Laboratory/NOAA, Princeton University,
!  Princeton N.J., U.S.A..
!
!  Killworth, P.D., Stainforth, D., Webb, D.J. and Paterson, P.M.,
!  1989:  A free surface Bryan-Cox-Semtner model.  Report No. 270,
!  Institute of Oceanographic Sciences, Wormley, U.K..
!
!  Pacanowski, R.C., Dixon, K., Rosati, A., 1990:  The GFDL Modular
!  Ocean Model 1.0.  Geophysical Fluid Dynamics Laboratory/NOAA,
!  Princeton University, Princeton, N.J., U.S.A..  (Unpublished
!  manuscript).
!
!  Webb, D.J., 1993:  An ocean model code for array processor
!  computers.  Internal Document No.324, Institute of Oceanographic
!  Sciences, Wormley, U.K..
!
!  Webb, D.J., 1993:  An ocean model code for array processor
!  computers.  Computers and Geosciences, 22(5), 569-578.
!
!==================================================================
!
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
!==================================================================
!
      real*4 utime(3)
!
      namelist /contrl/ init, fnrest, fnkmt, days, restrt, nmix, eb
     &,                 ncon, tsi, dgnstc, snapd, archd, acor
# ifdef de_checkbd
     &,                 dchkbd
# endif
     &,                 ftrest, ftsnap, ftarch
      namelist /tsteps/ dtts, dtuv, dtbt
      namelist /eddy/   am, ah, fkpm, fkph, cdbot
!
!  start timing
!
      call dtime(utime)
      utime(3)=0.0
!
!     constants usually initialised in block data
!
      model(1) =' Stripped down MOM code for array processor.'
      model(2) =' Module: moma.F.  Version: 1.17. '//
     &          ' Date: 09 Aug 1996'
      stlon = -4.0
      stlat = -72.0
      dxdeg = 4.0
      dydeg = 4.0
!
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
!
      init = .true.
      first= .true.
      days = 2.0
      eb    = .true.
      acor  = 0.6
      nmix  = 16
      ncon  = 1
      tsi   = 2.0
      dgnstc= 2.0
      snapd = 2.0
      archd = 2.0
      restrt = .false.
      fnrest = ' '
      fnkmt  = 'ocean.kmt'
      ftrest = 'hdf'
      ftsnap = 'hdf'
      ftarch = 'hdf'
      iorest = 21
!
      dtts = 10800d0
      dtuv = 10800d0
      dtbt = 100d0
!
      am    = 1.0e9
      ah    = 2.0e7
      fkph  = 20.0
      fkpm  = 1.0
      cdbot = 0.0
#ifdef de_checkbd
      dchkbd = 0.1
#endif
!
!==================================================================
!     begin introductory section
!     open input file and write MOM version information
!==================================================================
!
      open (stdin,file='ocean.in')
      write (stdout,'(/2(25x,a80/)/)') model
!
!------------------------------------------------------------------
!     initialize various quantities and
!     read in run parameters
!------------------------------------------------------------------
!
      pi     = c4*atan(c1)
      radian = c360/(c2*pi)
      omega  = pi/43082.0
      grav =980.6
      radius =6370.e5
!
      read  (stdin, contrl)
      write (stdout,contrl)
      read  (stdin, eddy)
      write (stdout,eddy)
      read  (stdin, tsteps)
      write (stdout,tsteps)
      close (unit=stdin,status='keep')
!
!------------------------------------------------------------------
!     calculate free surface model (barotropic) timestep
!------------------------------------------------------------------
!
      ntbt = max(1,nint(dtuv/dtbt))
      if(dtuv.ne.dtbt*ntbt)then
        write(stdout,'(a)') ' dtuv not integer multiple of dtbt.'
        write(stdout,*) ' dtuv is  ',dtuv
        write(stdout,*) ' dtbt was ',dtbt
        dtbt = dtuv/ntbt
        write(stdout,*) ' dtbt now set to set to:',dtbt
      endif
!
!------------------------------------------------------------------
!     list code options chosen by user for this model run
!------------------------------------------------------------------
!
      call docmnt (stdout, .false., .false., .false., .false.)
!
!------------------------------------------------------------------
!     set up model grids in x, y, & z
!------------------------------------------------------------------
!
      call grids
#ifdef de_checkbd
!
!------------------------------------------------------------------
! scale the checkerboard filter weighting factor (zw from grids)
!------------------------------------------------------------------
!
      dchkbd = dchkbd*grav*dtbt*zw(KM)
!
#endif
!
!==================================================================
!     end introductory section
!==================================================================
!
!------------------------------------------------------------------
!     set up initial pointers to timestep storage.  nc is current
!     timestep, np is next and nm is previous timestep.  nm0, nc0
!     and np0 are the corresponding pointers for the free surface
!     model.
!------------------------------------------------------------------
!
      nm = 1
      nc = 2
      np = 3
      nm0 = 1
      nc0 = 2
      np0 = 3
!
!------------------------------------------------------------------
!     is this a start from initial conditions or a restart?
!     in either case, set the number of vertical levels on the "t"
!     grid
!------------------------------------------------------------------
!
      if (init) then
        call ocn1st
      else
        call archrd
      endif
!
!------------------------------------------------------------------
!     compute number of vertical levels on the "u" grid
!------------------------------------------------------------------
!
      do 900 i=1,imt
        kmu(i,jmt) = 0
900   continue
      do 1000 j=1,jmtm1
      do 990 i=1,imtm1
        kmu(i,j) =min (kmt(i,j),kmt(i+1,j),kmt(i,j+1),kmt(i+1,j+1))
990   continue
1000  continue
      do 1100 j=1,jmt
#ifdef cyclic
        kmu(imt,j) = kmu(2,j)
#else
        kmu(imt,j) = 0
#endif
1100  continue
!
!------------------------------------------------------------------
!     compute area and volume of ocean ("t,s" grid boxes)
!------------------------------------------------------------------
!
      area   = c0
      volume = c0
      ocnp   = 0
      do 300 j=2,jmtm1
      do 300 i=2,imtm1
        if (kmt(i,j) .gt. 0) then
          area   = area + cst(j)*dx*dy
          volume = volume + cst(j)*dx*dy*zw(kmt(i,j))
          ocnp   = ocnp + float(kmt(i,j))
        endif
300   continue
      write (stdout,301) area, volume
!
!------------------------------------------------------------------
!       print map of "kmt" levels.
!------------------------------------------------------------------
!
!
      write (stdout,401)
      do 410 ibk=1,imt,40
        isp  = ibk
        iept = ibk + 40 - 1
        if(iept.gt.imt) iept=imt
        write (stdout,'(/, 4x, 40i3)') (ii, ii=isp,iept)
        do 400 jrev=1,jmt
          j=jmt-jrev+1
          write (stdout,'(1x,i3, 40i3)')j,(kmt(i,j),i=isp,iept)
 400    continue
 410  continue
      write(stdout,*)
!
!------------------------------------------------------------------
!     compute depths and reciprocal depths
!------------------------------------------------------------------
!
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
!
!------------------------------------------------------------------
!     initialize various things
!------------------------------------------------------------------
!
      do 2600 j=1,jmt
      do 2590 i=1,imt
        zu(i,j)  = c0
        zv(i,j)  = c0
2590  continue
2600  continue
!
!     initialize the coriolis factor
!
      do 2660 j=1,jmt
        fcor(j) = c2*omega*sine(j)
2660  continue
!
!     initialize the diffusion factors
!
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
!
!------------------------------------------------------------------
!     do all consistency checks now
!------------------------------------------------------------------
!
      call checks
!
!------------------------------------------------------------------
!     start the time step loop
!------------------------------------------------------------------
!
      do 3400 loop=1,9999999
        first= loop .eq. 1
!
!------------------------------------------------------------------
!       update timestep, set time dependent logical switches
!       to determine program flow for timestep itt, itt-1 and itt-2.
!       note: timestep itt and time refer to the time of the new
!       fields being calculated!
!------------------------------------------------------------------
!
        call tmngr (dtts)
        call step
!
!-----------------------------------------------------------------
!      print time step information, print diagnostics and archive
!      as required
!-----------------------------------------------------------------
!
        call dtime(utime)
        utime(3)=utime(1)+utime(2)+utime(3)
        if (prntsi) then
          write (stdout,9602) itt, stamp, ektot,
     &                      dtabs(1), dtabs(2),utime(1),utime(3)
          call flush(stdout)
        endif
        if (diagts) call pdiag
        if (snapts) call snap
        if (archts) call archwr
!
!------------------------------------------------------------------
!       end timestepping loop
!------------------------------------------------------------------
!
      if (last) go to 3401
3400  continue
3401  continue
!
!------------------------------------------------------------------
!     documentation info at end of run
!------------------------------------------------------------------
!
      call docmnt (stdout, .true., .true., .true., .true.)
!
!------------------------------------------------------------------
!     save restart dataset at end of run if restrt is true
!------------------------------------------------------------------
!
      if (restrt) call archwr
!
!------------------------------------------------------------------
!     close all units  (stderr is not closed because it is set equal
!     to stdout in pconst.h).
!------------------------------------------------------------------
!
      write (stdout,9603)
      close (unit=stdout,status='keep')
      stop
!
 301  format (//,'  Regional & Global ocean statistics:'
     &,/,'  the total ocean surface area (t grid) =',1pe15.8,'cm**2'
     &,/,'  the total ocean volume (t grid)       =',1pe15.8,'cm**3')
 401  format(/t50,'number of levels on "t,s" grid')
9451  format (/' ==== start and end indices for',a17,'====')
9461  format (' j=',i3,5x,5(2i5,10x))
9499  format (/' error => lseg too small for',a15,' indices'
     &        /'          j =',i5,'  lseg + 1 =',i8)
9602  format (1x,'ts=',i7, 1x, a32, ', ke=', 1pe13.6,
     &        ' dtemp=',1pe13.6,' dsalt=',1pe13.6,' times=',0p2f10.2)
9603  format(1x,' ==> END of model run.',//)
      end


      subroutine step
!
!==================================================================
!
!     step is called once per timestep. it includes all the main
!     loops over ic and jc and calls to the main routines.
!
!==================================================================
!
!
#include "param.h"
!
#include "scalar.h"
#include "switch.h"
#include "ctmngr.h"
#include "timelv.h"
#include "slabs.h"
#include "frees.h"
!
#include "cdiag.h"
!
!------------------------------------------------------------------
!    update pointers for new value of itt.
!    nnp, nnc and nnm are not changed during a timestep
!    np, nc, nm may be modified during a forward or euler backward
!    timestep.
!------------------------------------------------------------------
!
      nnc = np
      nnm = nc
      nnp = nm
      np  = nnp
      nc  = nnc
      nm  = nnm
!
!------------------------------------------------------------------
!       adjust various quantities for normal/mixing timesteps
!------------------------------------------------------------------
!
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
!
!------------------------------------------------------------------
!       return here for second pass of euler backward timestep
!------------------------------------------------------------------
!
 100  continue
#ifdef presetp
!
!-----------------------------------------------------------------
!
!     precalculate the baroclinic part of the pressure field
!     for use by subroutine clinic
!
!-----------------------------------------------------------------
!
      do 150 jc=1,jmt
      do 150 ic=1,imt
        call setp (ic,jc)
 150  continue
#endif
!
!-----------------------------------------------------------------
!      main baroclinic timestep loop over grid cells, from south
!      to north and from west to east
!      1. set vertical boundary conditions (surface & bottom)
!      2. calculate internal mode velocities
!      3. calculate tracers
!-----------------------------------------------------------------
!
      do 200 jc=2,jmtm1
      do 200 ic=2,imtm1
        call setvbc (ic,jc)
        call clinic (ic,jc)
        call tracer (ic,jc)
 200  continue
!
!-----------------------------------------------------------------
!      run free surface model (except during the second part
!      of a baroclinic euler backwards timestep).
!      first initialise pointers
!-----------------------------------------------------------------
!
      if(.not.mxpas2)then
        ntbtp0 = 0
        c2dtbt = dtbt
#ifdef free_eb
        ntbt2 = ntbt
#else
        ntbt2 = 2*ntbt
        do 230 jc=2,jmtm1
        do 230 ic=2,imtm1
        do 230 m=1,3
          freeav(m,ic,jc)=c0
 230    continue
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
!
!-----------------------------------------------------------------
!      On the first barotropic timestep or if "free_eb" is set
!      use euler backward. this requires two passes.
!-----------------------------------------------------------------
!
          do 500 itbtp=ntbtp0,1
!
!-----------------------------------------------------------------
!      main free surface model loop to carry out a
!      partial timestep for each model point
!-----------------------------------------------------------------
!
          do 300 jc=2,jmtm1
          do 300 ic=2,imtm1
            call frees(ic,jc)
 300      continue
!
!-----------------------------------------------------------------
!  set boundary conditions for the free surface model
!-----------------------------------------------------------------
!
          do 400 jc=2,jmtm1
#ifdef cyclic
            h0(  1,jc,np0) = h0(imtm1,jc,np0)
            h0(imt,jc,np0) = h0(    2,jc,np0)
            u0(  1,jc,np0) = u0(imum1,jc,np0)
            u0(imu,jc,np0) = u0(    2,jc,np0)
            v0(  1,jc,np0) = v0(imum1,jc,np0)
            v0(imu,jc,np0) = v0(    2,jc,np0)
#else
            h0(  1,jc,np0) = c0
            h0(imt,jc,np0) = c0
            u0(  1,jc,np0) = c0
            u0(imu,jc,np0) = c0
            v0(  1,jc,np0) = c0
            v0(imu,jc,np0) = c0
#endif
400       continue
!
!-----------------------------------------------------------------
!     reset pointers at end of first pass of the free surface model
!------------------------------------------------------------------
!
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
 500      continue
 600    continue
      endif
!
!-----------------------------------------------------------------
!  end of free surface model
!  now add barotropic velocities to baroclinic velocities
!-----------------------------------------------------------------
!
      do 700 jc=2,jmtm1
      do 700 ic=2,imtm1
        call addv(ic,jc)
 700  continue
!
!-----------------------------------------------------------------
!  set boundary conditions for the baroclinic model
!-----------------------------------------------------------------
!
      do 800 jc=2,jmtm1
      do 800 k =1,km
#ifdef cyclic
        u(k,  1,jc,np) = u(k,imum1,jc,np)
        u(k,imu,jc,np) = u(k,    2,jc,np)
        v(k,  1,jc,np) = v(k,imum1,jc,np)
        v(k,imu,jc,np) = v(k,    2,jc,np)
#else
        u(k,  1,jc,np) = c0
        u(k,imu,jc,np) = c0
        v(k,  1,jc,np) = c0
        v(k,imu,jc,np) = c0
#endif
        do 800 n=1,nt
#ifdef cyclic
          t(k,  1,jc,n,np) = t(k,imtm1,jc,n,np)
          t(k,imt,jc,n,np) = t(k,    2,jc,n,np)
#else
          t(k,  1,jc,n,np) = t(k,    1,jc,n,nm)
          t(k,imt,jc,n,np) = t(k,  imt,jc,n,nm)
#endif
 800  continue
!
!-----------------------------------------------------------------
!     if this is the end of the first pass of an euler backward
!     timestep then set the pointers for the second pass.
!-----------------------------------------------------------------
!
      if (mixts.and.eb)then
        eots   = .true.
        nc = nnp
        np = nnm
        mixts  = .false.
        mxpas2 = .true.
        go to 100
      endif
!
!-----------------------------------------------------------------
!  collect timestep statistics
!-----------------------------------------------------------------
!
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
!
!------------------------------------------------------------------
!       if this is the end of either a forward or euler backward
!       timestep then reset the pointers
!------------------------------------------------------------------
!
      if(mxpas2)then
        nc = nnc
        nm = nnp
      endif
      if(mixts)then
        nm = nnm
      endif
!
      return
      end


      subroutine clinic(ic,jc)
!
!==================================================================
!
!     clinic computes internal mode velocities at time level t+1
!     for a vertical column of k cells with i and j values equal to
!     ic and jc.
!
!==================================================================
!
!
#include "param.h"
!
#include "scalar.h"
#include "switch.h"
#include "coord.h"
#include "grdvar.h"
#include "levind.h"
#include "timelv.h"
#include "slabs.h"
#include "frees.h"
!
#include "chmix.h"
#include "cvbc.h"
#include "cvmix.h"
!
      kmc=kmu(ic,jc)
      if(kmc.eq.0)return
      jp = jc+1
      jm = jc-1
      ip = ic+1
      im = ic-1
!
!  find max of all the neighbouring points.
!
#ifdef oldadv
      kmd  = max(kmu(ip,jc),kmu(im,jc),kmu(ic,jp),kmu(ic,jm),kmc)
#else
      kmd  = max(kmu(im,jm),kmu(ic,jm),kmu(ip,jm),
     &           kmu(im,jc),kmc       ,kmu(ip,jc),
     &           kmu(im,jp),kmu(ic,jp),kmu(ip,jp))
#endif
      boxa = dx*csu(jc)*dy
      boxar = c1/boxa
!
!=============================================================
!  calculate horizontal and vertical advective fluxes
!=============================================================
!
      dxp5 = p5*dx
      dyp5 = p5*dy
      do 100 k=1,kmd
#ifdef oldadv
!
!  original horizontall advection scheme
!
      fvn(k) = (v(k,ic,jc,nc)+v(k,ic,jp,nc))*p5*dx*cst(jp)*dz(k)
      fvs(k) = (v(k,ic,jc,nc)+v(k,ic,jm,nc))*p5*dx*cst(jc)*dz(k)
      fue(k) = (u(k,ic,jc,nc)+u(k,ip,jc,nc))*p5*dy*dz(k)
      fuw(k) = (u(k,ic,jc,nc)+u(k,im,jc,nc))*p5*dy*dz(k)
#else
!
!  new horizontal advection scheme
!
      fvjc   = csu(jc)*
     &         (v(k,im,jc,nc)+2.0*v(k,ic,jc,nc)+v(k,ip,jc,nc))
      fuic   = (u(k,ic,jp,nc)+2.0*u(k,ic,jc,nc)+u(k,ic,jm,nc))
!
      fvn(k) = dz(k)*dx*p125*( fvjc + csu(jp)*
     &         (v(k,im,jp,nc)+2.0*v(k,ic,jp,nc)+v(k,ip,jp,nc)))
      fvs(k) = dz(k)*dx*p125*( fvjc + csu(jm)*
     &         (v(k,im,jm,nc)+2.0*v(k,ic,jm,nc)+v(k,ip,jm,nc)))
!
      fue(k) = dz(k)*dy*p125*( fuic +
     &         (u(k,ip,jp,nc)+2.0*u(k,ip,jc,nc)+u(k,ip,jm,nc)))
      fuw(k) = dz(k)*dy*p125*( fuic +
     &         (u(k,im,jp,nc)+2.0*u(k,im,jc,nc)+u(k,im,jm,nc)))
#endif
!
!  new vertical boundary conditions on w at u points
!
      fw(k-1)  = (fuw(k)-fue(k)+fvs(k)-fvn(k))
 100  continue
!
!=============================================================
!    integrate up from the bottom level
!=============================================================
!
      fw(kmd)=c0
      do 210 k=kmd,1,-1
      fw(k-1)=fw(k-1)+fw(k)
 210  continue
#ifdef presetp
!
!-----------------------------------------------------------------
!         compute hydrostatic pressure gradient (baroclinic part)
!-----------------------------------------------------------------
!
      fxa = dx2r*csur(jc)
      fxb = dy2r
      do 300 k=1,kmc
        temp1  = p(k,ip,jp) - p(k,ic,jc)
        temp2  = p(k,ic,jp) - p(k,ip,jc)
        dpdx(k) = (temp1-temp2)*fxa
        dpdy(k) = (temp1+temp2)*fxb
 300  continue
#else
!
!-----------------------------------------------------------------
!     compute density on surrounding 't' points
!-----------------------------------------------------------------
!
      call state (t(1,ic,jc,1,nc), t(1,ic,jc,2,nc), rhoo,
     &            kmt(ic,jc) ,0)
      call state (t(1,ip,jc,1,nc), t(1,ip,jc,2,nc), rhpo,
     &            kmt(ip,jc) ,0)
      call state (t(1,ip,jp,1,nc), t(1,ip,jp,2,nc), rhpp,
     &            kmt(ip,jp) ,0)
      call state (t(1,ic,jp,1,nc), t(1,ic,jp,2,nc), rhop,
     &            kmt(ic,jp) ,0)
!
!-----------------------------------------------------------------
!         compute hydrostatic pressure gradient (baroclinic part)
!-----------------------------------------------------------------
!
!    1. compute it at the first level
!
      fxa = grav*dzw(0)*dx2r*csur(jc)
      fxb = grav*dzw(0)*dy2r
      temp1  = rhpp(1) - rhoo(1)
      temp2  = rhop(1) - rhpo(1)
      dpdx(1) = (temp1-temp2)*fxa
      dpdy(1) = (temp1+temp2)*fxb
!
!    2. compute the change in pressure gradient between levels
!
      fxa = grav*dx4r*csur(jc)
      fxb = grav*dy4r
      do 300 k=2,kmc
        temp1 = rhpp(k)+rhpp(k-1)- (rhoo(k)+rhoo(k-1))
        temp2 = rhop(k)+rhop(k-1)- (rhpo(k)+rhpo(k-1))
        dpdx(k) = (temp1-temp2)*fxa*dzw(k-1)
        dpdy(k) = (temp1+temp2)*fxb*dzw(k-1)
 300  continue
!
!=============================================================
!         integrate downward from the first level
!=============================================================
!
      do 400 k=2,kmc
        dpdx(k) = dpdx(k-1) + dpdx(k)
        dpdy(k) = dpdy(k-1) + dpdy(k)
 400  continue
#endif
!
!-----------------------------------------------------------------
!         calculate quantities for the computation of
!         vertical diffusion and advection of momentum
!         - positive for upwards flux of momentum
!-----------------------------------------------------------------
!
      do 500 k=1,kmc-1
        vmf(k,1) = fkpm*(u(k+1,ic,jc,nm) - u(k  ,ic,jc,nm))
     &                     *dzwr(k)*dx*csu(jc)*dy
        vmf(k,2) = fkpm*(v(k+1,ic,jc,nm) - v(k  ,ic,jc,nm))
     &                     *dzwr(k)*dx*csu(jc)*dy
        fwb1(k) = fw(k)*(u(k,ic,jc,nc) + u(k+1,ic,jc,nc))*p5
        fwb2(k) = fw(k)*(v(k,ic,jc,nc) + v(k+1,ic,jc,nc))*p5
 500  continue
!
!-----------------------------------------------------------------
!     set k=0 vertical diffusion to reflect surface wind stress
!     set level "kmc" for bottom drag condition.
!     set k=kmc elements of vertical advection array to zero
!     vertical momentum transfer.  Set k=0 to represent advection
!     loss to region above the free surface.
!     smf is the surface stress acting on the ocean
!     bmf is the bottom stress acting on the ocean floor
!-----------------------------------------------------------------
!
      vmf(0,1)  = -smf(1)*boxa
      vmf(0,2)  = -smf(2)*boxa
      vmf(kmc,1)= -bmf(1)*boxa
      vmf(kmc,2)= -bmf(2)*boxa
      fwb1(0)   = fw(0)*u(1,ic,jc,nc)
      fwb2(0)   = fw(0)*v(1,ic,jc,nc)
      fwb1(kmc) = c0
      fwb2(kmc) = c0
!
!------------------------------------------------------------------
!       set up index lcor for the coriolis term:
!       use tau timestep for explicit trtmnt ... or
!       use tau-1 timestep for implicit treatment with remainder of
!       term to be added later
!------------------------------------------------------------------
!
          if (acor .eq. c0) then
            lcor = nc
          else
            lcor = nm
          endif
!
!------------------------------------------------------------------
!       initialise vertical mean forcing array elements:
!------------------------------------------------------------------
!
      zu(ic,jc) = c0
      zv(ic,jc) = c0
!
!=============================================================
!    main calculation to timestep baroclinic velocity
!=============================================================
!
      do 600 k=1,kmc
!
! u-momentum equation.  separate coriolis term and the rest.
!
      ucor  = fcor(jc)*v(k,ic,jc,lcor)
      urest = -dpdx(k)
!  u-advection
     & +(fuw(k)*(u(k,im,jc,nc) + u(k,ic,jc,nc))*p5
     & - fue(k)*(u(k,ip,jc,nc) + u(k,ic,jc,nc))*p5
!  v-advection
     & + fvs(k)*(u(k,ic,jm,nc) + u(k,ic,jc,nc))*p5
     & - fvn(k)*(u(k,ic,jp,nc) + u(k,ic,jc,nc))*p5
!  w-advection
     & + fwb1(k  ) - fwb1(k-1)
!  x diffusion
     & +(bbu(jc)*(u(k,ip,jc,nm)-u(k,ic,jc,nm))
     & + bbu(jc)*(u(k,im,jc,nm)-u(k,ic,jc,nm))
!  y diffusion
     & + ccu(jc)*(u(k,ic,jp,nm)-u(k,ic,jc,nm))
     & + ddu(jc)*(u(k,ic,jm,nm)-u(k,ic,jc,nm)) )
     & *dz(k)
!  z diffusion
     &  + vmf(k  ,1) - vmf(k-1,1) )*boxar*dzr(k)
!
! v-momentum equation.  separate coriolis term and the rest
!
      vcor  = -fcor(jc)*u(k,ic,jc,lcor)
      vrest = -dpdy(k)
!  u-advection
     & +(fuw(k)*(v(k,im,jc,nc) + v(k,ic,jc,nc))*p5
     & - fue(k)*(v(k,ip,jc,nc) + v(k,ic,jc,nc))*p5
!  v-advection
     & + fvs(k)*(v(k,ic,jm,nc) + v(k,ic,jc,nc))*p5
     & - fvn(k)*(v(k,ic,jp,nc) + v(k,ic,jc,nc))*p5
!  w-advection
     & + fwb2(k  ) - fwb2(k-1)
!  x diffusion
     & +(bbu(jc)*(v(k,ip,jc,nm)-v(k,ic,jc,nm))
     & + bbu(jc)*(v(k,im,jc,nm)-v(k,ic,jc,nm))
!  y diffusion
     & + ccu(jc)*(v(k,ic,jp,nm)-v(k,ic,jc,nm))
     & + ddu(jc)*(v(k,ic,jm,nm)-v(k,ic,jc,nm)) )
     & *dz(k)
!  z diffusion
     & + vmf(k  ,2) - vmf(k-1,2) )*boxar*dzr(k)
!
!***************************************************************
!  the corrections to the advection and diffusion terms that
!  allow for the curvature of the grid near the poles have been
!  left out of this version of the model
!***************************************************************
!
!-----------------------------------------------------------------
!     calculate vertically averaged forcing (without coriolis)
!     and total force at each level (with coriolis)
!-----------------------------------------------------------------
!
      fx = dz(k)
      zu(ic,jc) = zu(ic,jc) + urest*fx
      zv(ic,jc) = zv(ic,jc) + vrest*fx
      u(k,ic,jc,np) = ucor + urest
      v(k,ic,jc,np) = vcor + vrest
 600  continue
!
!-----------------------------------------------------------------
!     calculate vertical mean force
!-----------------------------------------------------------------
!
      zu(ic,jc) = zu(ic,jc)*hr(ic,jc)
      zv(ic,jc) = zv(ic,jc)*hr(ic,jc)
!
!-----------------------------------------------------------------
!     timestep baroclinic terms
!-----------------------------------------------------------------
!
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
!
!-----------------------------------------------------------------
!     remove vertical mean velocity
!-----------------------------------------------------------------
!
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
!
!=============================================================
!
!     tracer computes tracer quantities for time level t+1
!
!=============================================================
!
#include "param.h"
!
#include "scalar.h"
#include "grdvar.h"
#include "levind.h"
#include "timelv.h"
#include "slabs.h"
!
#include "chmix.h"
#include "cvbc.h"
#include "cvmix.h"
!
!=============================================================
!    outer loops over i and j indices
!=============================================================
!
      kmc= kmt(ic,jc)
      if(kmc.eq.0)return
      jp = jc+1
      jm = jc-1
      ip = ic+1
      im = ic-1
      boxa = dx*cst(jc)*dy
      boxar = c1/boxa
!
!=============================================================
!    set up masks, equal to 0 on sub-seabottom points.
!=============================================================
!
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
!
!=============================================================
!  calculate horizontal and vertical advective fluxes
!=============================================================
!
      dxp5 = p5*dx
      dyp5 = p5*dy
      do 160 k=1,kmc
      fvn(k) = (v(k,ic,jc,nc)+v(k,im,jc,nc))*dxp5*csu(jc)*dz(k)
      fvs(k) = (v(k,ic,jm,nc)+v(k,im,jm,nc))*dxp5*csu(jm)*dz(k)
      fue(k) = (u(k,ic,jc,nc)+u(k,ic,jm,nc))*dyp5*dz(k)
      fuw(k) = (u(k,im,jc,nc)+u(k,im,jm,nc))*dyp5*dz(k)
      fw(k-1)  = (fuw(k)-fue(k)+fvs(k)-fvn(k))
 160  continue

!
!=============================================================
!    integrate up from the lowest layer
!    fw is positive upwards
!=============================================================
!
      fw(kmc)=c0
      do 180 k=kmc-1,0,-1
      fw(k)=fw(k)+fw(k+1)
 180  continue
!
!=============================================================
!    main calculation loop to timestep tracers
!=============================================================
!
      do 900 m=1,nt
!
!=============================================================
!  calculate vertical fluxes - positive upwards
!=============================================================
!
      do 930 k=1,kmc-1
        vtf(k) = fkph*(t(k+1,ic,jc,m,nm)-t(k  ,ic,jc,m,nm))
     &                 *dzwr(k)*boxa
        fwb1(k) = fw(k)*(t(k,ic,jc,m,nc)+t(k+1,ic,jc,m,nc))*p5
 930  continue
!
!=============================================================
!  surface and bottom advective and diffusive boundary conditions
!  stf is downwards flux into the ocean
!=============================================================
!
      vtf(0)   = -stf(m)*boxa
      vtf(kmc) = c0
      fwb1(0)  = fw(0)*t(1,ic,jc,m,nc)
      fwb1(kmc)= c0
!
!=============================================================
!  main timestepping calculation
!=============================================================
!
!
      do 960 k=1,kmc
      t(k,ic,jc,m,np) = t(k,ic,jc,m,nm)
     &                + c2dtts*boxar*dzr(k)*(
!  u-advection
     &   fuw(k)*(t(k,im,jc,m,nc) + t(k,ic,jc,m,nc))*p5
     & - fue(k)*(t(k,ip,jc,m,nc) + t(k,ic,jc,m,nc))*p5
!  v-advection
     & + fvs(k)*(t(k,ic,jm,m,nc) + t(k,ic,jc,m,nc))*p5
     & - fvn(k)*(t(k,ic,jp,m,nc) + t(k,ic,jc,m,nc))*p5
!  w-advection
     & + fwb1(k  ) - fwb1(k-1)
!  x diffusion
     & +(bbt(jc)*(t(k,ip,jc,m,nm)-t(k,ic,jc,m,nm))*maskpo(k)
     & + bbt(jc)*(t(k,im,jc,m,nm)-t(k,ic,jc,m,nm))*maskmo(k)
!  y diffusion
     & + cct(jc)*(t(k,ic,jp,m,nm)-t(k,ic,jc,m,nm))*maskop(k)
     & + ddt(jc)*(t(k,ic,jm,m,nm)-t(k,ic,jc,m,nm))*maskom(k))
     & * dz(k)
!  z diffusion
     & + vtf(k  ) - vtf(k-1)
     &   )
 960  continue
 900  continue
!
!=============================================================
!   convectively adjust water column if gravitaionally unstable.
!=============================================================
!
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
!
      return
      end


      subroutine state (t, s, rho, kk, ind)
!
!=============================================================
!
!     state computes one column of normalized densities from the
!     surface to level kk by using a 3rd order polynomial fit to
!     the equation of state.
!
!     note.. for precision purposes, there is a depth dependent
!     constant subtracted from the density returned by this
!     routine.
!     so.. this routine should be used only for calculating
!     horizontal gradients of density or vertical stability
!     with adjacent levels referred to the same reference level.
!
!     inputs:
!
!     t = the input column of temperatures (degrees C)
!     s = the input column of salinities (units: (ppt-35)/1000)
!     ind = 0 for normal density calculations
!         = 1 for comparing levels 1 to 2, 3 to 4, etc.
!         = 2 for comparing levels 2 to 3, 4 to 5, etc.
!           (if ind = 1 or 2, the coefficients for the lower
!            of the 2 levels are used)
!     kk = bottom level to be calculated
!
!     note: model salinities are equal to (S-35.0)/1000.0 where
!     S is the standard oceanographic salinity.  This convention,
!     which dates from the early Bryan and Cox model, increases
!     the precision by about 5 binary bits.  This can be
!     significant in density gradient calculations on 32-bit
!     computers.  The convention is also the one used by the
!     mom/gfdl program used to generate the include file
!     'dncoef.h'.
!
!     output:
!
!     rho = density with depth level dependent offset
!           (units: g/cc minus constant)
!
!=============================================================
!
!
#include "param.h"
      dimension t(km),  s(km),  rho(km)
      dimension to(km), so(km), c(km,9),
     &          to1(km),so1(km),co1(km,9),
     &          to2(km),so2(km),co2(km,9)
      save in,to,so,c,to1,so1,co1,to2,so2,co2
#include "dncoef.h"
      data in/1/
!
!  on first entry set up coefficients for stability calculations
!
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
!
      do 210 k=1,km
        rho(k)=1e30
 210  continue
!
!  test for errors
!
      if(ind.lt.0.or.ind.gt.2.or.kk.lt.0.or.kk.gt.km)then
        write (stderr,99)ind,kk
        stop '1 state'
!
!  normal entry
!
      elseif(ind.eq.0)then
        do 400 k=1,kk
        tq=t(k)-to(k)
        sq=s(k)-so(k)
        rho(k) = (c(k,1) + (c(k,4) + c(k,7)*sq)*sq +
     &           (c(k,3) + c(k,8)*sq + c(k,6)*tq)*tq)*tq +
     &           (c(k,2) + (c(k,5) + c(k,9)*sq)*sq)*sq

 400    continue
!
!  compare 1 with 2, 3 with 4, etc.  density w.r.t bottom level
!
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
!
!  compare 2 with 3, 4 with 5 etc
!
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
!
      return
   99 format(/,' Routine state. Variables "ind" or "kk" out',
     &       ' of range',/,' ind =',i10,' kk =',i10)
      end



      subroutine frees(ic,jc)
!
!====================================================================
!
!  Based on the routine in "Killworth, P.D., Stainforth, D.,
!  Webb, D.J. and Paterson, S.M. (1989).  A free surface
!  Bryan-Cox-Semtner model.  Institute of Oceanographic
!  Sciences, Report No. 270.  Wormley, Godlaming, U.K.. 184pp."
!
!  Later modified to use a leapfrog and time smoothing scheme,
!  suggested by Semtner, to prevent aliasing instabilities
!  near the equator
!
!=====================================================================
!
!---------------------------------------------------------------------
!  define global data
!---------------------------------------------------------------------
!
#include "param.h"
!
#include "chmix.h"
#include "ctmngr.h"
#include "frees.h"
#include "grdvar.h"
#include "levind.h"
#include "scalar.h"
#include "switch.h"
#include "timelv.h"
!
!
!---------------------------------------------------------------------
!  set indices
!---------------------------------------------------------------------
!
      ip=ic+1
      im=ic-1
      jp=jc+1
      jm=jc-1
!
!---------------------------------------------------------------------
!  timestep height field.
!---------------------------------------------------------------------
!
      if(kmt(ic,jc).ne.0)then
        boxar = dxr*dyr*cstr(jc)
#ifdef de_checkbd
!
! Calculate masks for the del-cross and del-plus operators
! on selected time steps for free-surface calculation
!
       lchkbd = itbt.eq.2
       if (lchkbd) then
         m_N  = min(1,kmt(ic,jp))
         m_S  = min(1,kmt(ic,jm))
         m_E  = min(1,kmt(ip,jc))
         m_W  = min(1,kmt(im,jc))
         m_NE = min(1,kmu(ic,jc))
         m_SE = min(1,kmu(ic,jm))
         m_SW = min(1,kmu(im,jm))
         m_NW = min(1,kmu(im,jc))
!
! Calculate the del-plus operator (x2 for inclusion within existing terms
! of the dhdt expression)
!
         delplus = c2*(m_E*(h0(ip,jc,nm0) - h0(ic,jc,nm0))
     &                -m_W*(h0(ic,jc,nm0) - h0(im,jc,nm0))
     &                +m_N*(h0(ic,jp,nm0) - h0(ic,jc,nm0))
     &                -m_S*(h0(ic,jc,nm0) - h0(ic,jm,nm0)))
!
! Calculate the del-cross operator (x2 for inclusion within existing terms
! of the dhdt expression)
!
         delcross = ( m_NE*(h0(ip,jp,nm0) - h0(ic,jc,nm0))
     &               -m_SW*(h0(ic,jc,nm0) - h0(im,jm,nm0))
     &               +m_NW*(h0(im,jp,nm0) - h0(ic,jc,nm0))
     &               -m_SE*(h0(ic,jc,nm0) - h0(ip,jm,nm0)))
       endif
!
! Need a constant of the form:
!        WGHT = alpha*grav*Hmax*dtbt/(dphi*dlambda*radius*radius*cos(phi))
!    i.e.WGHT = alpha*grav*Hmax*dtbt*boxar = dchkbd*boxar
!
# endif
        dhdt  = (dy*(u0(im,jm,nc0)*h(im,jm)+u0(im,jc,nc0)*h(im,jc)
     &              -u0(ic,jm,nc0)*h(ic,jm)-u0(ic,jc,nc0)*h(ic,jc))
     & +dx*(csu(jm)*(v0(im,jm,nc0)*h(im,jm)+v0(ic,jm,nc0)*h(ic,jm))
     &     -csu(jc)*(v0(im,jc,nc0)*h(im,jc)+v0(ic,jc,nc0)*h(ic,jc))
     &     ))*0.5*boxar
# ifdef de_checkbd
       if(lchkbd) then
        dhdt=dhdt+dchkbd*(delplus-delcross)*0.5*boxar
       endif
# endif
!
        h0(ic,jc,np0) = h0(ic,jc,nm0) + c2dtbt*dhdt
#ifndef free_eb
!
!---------------------------------------------------------------------
!  Calculate time averages (except for first pass of first timestep
!                           when itbtp equals zero)
!---------------------------------------------------------------------
!
        if(itbtp.ne.0)then
          freeav(1,ic,jc) = freeav(1,ic,jc) + h0(ic,jc,np0)
          if(itbt.eq.ntbt2)then
            fx = 0.5/ntbt
            h0(ic,jc,np0) = freeav(1,ic,jc)*fx
          endif
        endif
#endif
      endif
!
!---------------------------------------------------------------------
!  timestep velocity fields using semi-implicit scheme
!  for the coriolis term.
!---------------------------------------------------------------------
!
      if(kmu(ic,jc).ne.0)then
!
!---------------------------------------------------------------------
!    compute pressure gradients
!---------------------------------------------------------------------
!
        fxa = grav*dx2r*csur(jc)
        fxb = grav*dy2r
        temp1 = h0(ip,jp,nc0) - h0(ic,jc,nc0)
        temp2 = h0(ic,jp,nc0) - h0(ip,jc,nc0)
        dpdx  = (temp1-temp2)*fxa
        dpdy  = (temp1+temp2)*fxb
#ifdef free_eb
!
!---------------------------------------------------------------------
!  euler backward scheme treats coriolis term implicitly
!---------------------------------------------------------------------
!
        fac1=fcor(jc)*c2dtbt*p5
        fac2=c1/(c1+fac1*fac1)
        ustar = u0(ic,jc,nm0) + c2dtbt*(zu(ic,jc) - dpdx)
     &                        + fac1*v0(ic,jc,nm0)
        vstar = v0(ic,jc,nm0) + c2dtbt*(zv(ic,jc) - dpdy)
     &                        - fac1*u0(ic,jc,nm0)
        u0(ic,jc,np0) = (ustar + fac1*vstar)*fac2
        v0(ic,jc,np0) = (vstar - fac1*ustar)*fac2
#else
!
!---------------------------------------------------------------------
!  leapfrog scheme (explicit matsuno/euler-backward
!                   on first mixing timestep)
!---------------------------------------------------------------------
!
        u0(ic,jc,np0) = u0(ic,jc,nm0) + c2dtbt*(zu(ic,jc) - dpdx
     &                        + fcor(jc)*v0(ic,jc,nc0) )
        v0(ic,jc,np0) = v0(ic,jc,nm0) + c2dtbt*(zv(ic,jc) - dpdy
     &                        - fcor(jc)*u0(ic,jc,nc0) )
!
!---------------------------------------------------------------------
!  Calculate time averages (except for first pass of matsuno/
!                           euler-backward when itbtp equals zero)
!---------------------------------------------------------------------
!
        if(itbtp.ne.0)then
          freeav(2,ic,jc) = freeav(2,ic,jc) + u0(ic,jc,np0)
          freeav(3,ic,jc) = freeav(3,ic,jc) + v0(ic,jc,np0)
          if(itbt.eq.ntbt2)then
            fx = 0.5/ntbt
            u0(ic,jc,np0) = freeav(2,ic,jc)*fx
            v0(ic,jc,np0) = freeav(3,ic,jc)*fx
          endif
        endif
#endif
      endif
      return
      end


#ifdef presetp
      subroutine setp (ic, jc)
!
!==================================================================
!
!     precalculate the baroclinic part of the pressure field
!     for use by subroutine clinic
!
!     p(1,...) statement moved within loop 100 to overcome bug
!     with Silicon Graphics at level 3 optimisation
!
!==================================================================
!
#include "param.h"
#include "grdvar.h"
#include "scalar.h"
#include "slabs.h"
#include "timelv.h"
#include "levind.h"
!
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
      return
      end
#endif



      subroutine setvbc (ic, jc)
!
!==================================================================
!
!     compute vertical boundary conditions for clinic & tracer
!
!==================================================================
!
#include "param.h"
#include "cvbc.h"
#include "grdvar.h"
#include "scalar.h"
#include "slabs.h"
#include "timelv.h"
#include "levind.h"
!
!-----------------------------------------------------------------
!     get test case surface temperatures, salinities and wind
!     stresses (they vary with latitude).  test assumes nt=2,
!     (potential temperature and salinity)
!------------------------------------------------------------------
!
      call bcest(stf,smf, ic, jc)
!  convert salinity to model units.
      stf(2) = (stf(2)-35.0)*0.001
!
!------------------------------------------------------------------
!     set no flux condition for tracer at  bottom.
!     relax to "bcest" values with 50 day timescale at top
!------------------------------------------------------------------
!
      gamma = secday/50.0
      do 100 m=1,nt
        stf(m) = gamma*(stf(m)-t(1,ic,jc,m,nm))*dz(1)
        btf(m) = c0
100   continue
!
!------------------------------------------------------------------
!     set bottom stresses.  bmf is stress acting on the ocean floor
!------------------------------------------------------------------
!
      kz = kmu(ic,jc)
      if (kz .ne. 0) then
        uvmag    = sqrt(u(kz,ic,jc,nm)**2 +v(kz,ic,jc,nm)**2)
        bmf(1) = cdbot*u(kz,ic,jc,nm)*uvmag
        bmf(2) = cdbot*v(kz,ic,jc,nm)*uvmag
      else
        bmf(1) = c0
        bmf(2) = c0
      endif
!
      return
      end


      subroutine bcest (stf, smf, i, j)
!
!==================================================================
!     this routine linearly interpolates global, zonal mean values
!     of ocean model surface boundary conditions (SST, salinity,
!     WSX, WSY) having 4.5 deg spacing, to the requested latitdude.
!==================================================================
#include "param.h"
      parameter (nolat = 40, nolatp = nolat+1)
#include "coord.h"
      dimension sstobs(nolat), salobs(nolat)
      dimension wsxobs(nolatp), wsyobs(nolatp)
      dimension olatt(nolat), olatv(nolatp)
      dimension stf(nt), smf(2)
!
!     units of boundary condition values are:
!       degrees C, parts per thousand, dynes per square centimeter
!     olatt  = latitude points for observed data
!     olatv  = latitude points for observed data
!     dolat  = latitude spacing for observed data
!
!     "observed" temperature and salinity data are based on global,
!     annual mean zonally averaged values from the Levitus Atlas
!     (1982).
!     "observed" windstress data are based on global, annual mean,
!     zonally averaged values from Hellerman and Rosenstein (1981).
!     some smoothing was done.
!
!     units:  temperature - degrees Centigrade
!             salinity    - standard oceanographic salinity units
!             wind stress - dynes per sq cm
!     The salinity is transformed into model units in routine
!     'setvbc'.
!
!     references:
!       Hellerman, S, and M. Rosenstein, normal monthly wind stress
!         over the world ocean with error estimates, Journal of
!         Physical Oceanography, 13, 1093-1104,1983.
!       Levitus, S., Climatological atlas of the world ocean, NOAA
!         Prof. Paper 13, US Gov't printing Office, Washington, DC,
!         1982.
!
      save in,dolat,rdolat,olatt,olatv,sstobs,salobs,wsxobs,wsyobs
!
      data sstobs / -1.75, -1.75, -1.50, -1.50, -1.28,
     &              -0.55,  0.90,  2.92,  5.45,  8.62,
     &              12.27, 15.49, 18.30, 20.67, 22.64,
     &              24.14, 25.27, 26.37, 26.52, 26.16,
     &              26.85, 27.27, 26.82, 26.42, 25.53,
     &              24.03, 22.07, 19.73, 17.02, 12.77,
     &               8.93,  7.25,  6.22,  4.67,  4.57,
     &               3.03, -0.01, -1.05, -1.75, -1.75/
!
      data salobs / 34.30, 34.30, 34.30, 34.13, 33.98,
     &              33.97, 33.97, 33.98, 34.03, 34.24,
     &              34.61, 35.02, 35.37, 35.61, 35.72,
     &              35.68, 35.51, 35.22, 35.05, 35.12,
     &              34.80, 34.56, 34.71, 34.90, 35.27,
     &              35.67, 35.56, 35.49, 35.23, 34.28,
     &              33.57, 33.57, 33.60, 33.80, 34.04,
     &              34.05, 32.65, 32.30, 32.10, 32.00/
!
      data wsxobs /  0.00,
     &               0.00,  0.00, -0.02,  0.15,  0.31,
     &               0.50,  0.82,  1.08,  1.23,  1.16,
     &               0.84,  0.41,  0.02, -0.35, -0.55,
     &              -0.67, -0.64, -0.46, -0.29, -0.19,
     &              -0.16, -0.33, -0.52, -0.59, -0.55,
     &              -0.32,  0.09,  0.42,  0.56,  0.76,
     &               0.81,  0.65,  0.29,  0.06, -0.10,
     &              -0.05, -0.03,  0.05,  0.10,  0.01/
!
      data wsyobs /  .000,
     &               .000,  .009,  .032,  .005, -.023,
     &              -.075, -.155, -.202, -.230, -.179,
     &              -.049,  .093,  .214,  .294,  .344,
     &               .383,  .364,  .269,  .189,  .178,
     &               .125, -.122, -.213, -.251, -.259,
     &              -.202, -.189, -.179, -.183, -.009,
     &               .023,  .053, -.048, -.185, -.225,
     &              -.097, -.050, -.023, -.006,  .000/
!
      data in/1/
!
!----------------------------------------------------------------
!     set latitudes of sst and salinity observations
!     and set latitudes of windstress observations
!----------------------------------------------------------------
!
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
!
!----------------------------------------------------------------
!   use linear interpolation to produce the estimated surface
!   boundary condition values for temperature and salinity at t,s
!   row j
!----------------------------------------------------------------
!
!
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
!
!----------------------------------------------------------------
!   use linear interpolation to produce the estimated surface
!   boundary condition values for wind stress components at
!   u,v row j
!----------------------------------------------------------------
!
      ylatv = yu(j)
!
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
!
      return
      end


      subroutine addv( ic, jc)
!
!==================================================================
!
!  add barotropic to baroclinic velocities
!
!==================================================================
!
#include "param.h"
!
#include "frees.h"
#include "levind.h"
#include "timelv.h"
#include "slabs.h"
!
      ubar=u0(ic,jc,np0)
      vbar=v0(ic,jc,np0)
      do 100 k=1,kmu(ic,jc)
        u(k,ic,jc,np) = u(k,ic,jc,np) + ubar
        v(k,ic,jc,np) = v(k,ic,jc,np) + vbar
 100  continue
!
      return
      end


      subroutine diag( ic, jc)
!
!=============================================================
!  collect timestep statistics and (if diagts is set) statistics
!  for diagnostic timestep.
!  note 1. this routine sums arrays, it may not vectorise well
!       2. diagnostic timestep option not implemented yet
!=============================================================
!
!-----------------------------------------------------------------
!         accumulate kinetic energy
!-----------------------------------------------------------------
!
!
#include "param.h"
!
#include "scalar.h"
#include "switch.h"
#include "grdvar.h"
#include "levind.h"
#include "timelv.h"
#include "slabs.h"
!
#include "cdiag.h"
!
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
!
!-----------------------------------------------------------------
!         accumulate tracer change statistics
!-----------------------------------------------------------------
!
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
!==================================================================
!
!     calculate all grid related items
!
!==================================================================
!
#include "param.h"
#include "scalar.h"
#include "coord.h"
#include "grdvar.h"
!
!------------------------------------------------------------------
!     set latitudes & longitudes of "t" & "u,v" grid points (in degrees)
!------------------------------------------------------------------
!
      do 100 i=1,imt
        xt(i) = stlon + (i-p5)*dxdeg
        xu(i) = stlon + i*dxdeg
100   continue
!
      do 200 j=1,jmt
        yt(j) = stlat + (j-p5)*dydeg
        yu(j) = stlat + j*dydeg
200   continue
!
!     convert to centimeters
!
      dy = dydeg*radius/radian
      dx = dxdeg*radius/radian
!
!------------------------------------------------------------------
!     compute auxiliary arrays
!------------------------------------------------------------------
!
      do 500 k=1,km
        dzr(k)  = c1/dz(k)
        c2dz(k) = c2*dz(k)
        dz2r(k) = c1/c2dz(k)
500   continue
!
      dzw(0) = p5*dz(1)
      zw(1)  = dz(1)
      do 600 k=2,km
        dzw(k-1) = p5*(dz(k-1)+dz(k))
        zw(k)    = zw(k-1)+dz(k)
600   continue
!
      dzw(km)   = p5*dz(km)
      dzwr(km)  = c1/dzw(km)
      dzw2r(km) = p5/dzw(km)
      zt(1)     = dzw(0)
!
      do 700 k=1,km
        dzwr(k-1)  = c1/dzw(k-1)
        dzw2r(k-1) = p5/dzw(k-1)
        zt(k+1)  = zt(k) + dzw(k)
700   continue
!
      dxr = c1/dx
      dx2r = p5/dx
      dx4r = p25/dx
!
      dyr = c1/dy
      dy2r = p5/dy
      dy4r = p25/dy
!
      phi(1)   = yu(1)/radian
      phit(1)  = phi(1) - p5*dydeg/radian
      sumdy    = phi(1)
!
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
!
!------------------------------------------------------------------
!     print grid related arrays
!------------------------------------------------------------------
!
!
      write (stdout,9101)
      write (stdout,9001) dz
      write (stdout,9102)
      write (stdout,9001) dzw
      write (stdout,9103)
      write (stdout,9001) zw
      write (stdout,9104)
      write (stdout,9001) zt
      write (stdout,9105)
      write (stdout,9001) yt
      write (stdout,9106)
      write (stdout,9001) yu
      write (stdout,9107)
      write (stdout,9001) xt
      write (stdout,9108)
      write (stdout,9001) xu
!
      return
!
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
!
!==================================================================
!
!  initialize quantites when model first started
!
!==================================================================
!
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
!
!------------------------------------------------------------------
!     generate topography over "t" points
!------------------------------------------------------------------
!
#ifdef simple_world
      call topog
#else
      call read_kmt
#endif
!
!------------------------------------------------------------------
!     initialize t,s,u,v.
!          velocities set to zero
!          temperature an salinity set to 0.0 and 35.0 at sea
!                                      to -2.0 and 45.0 over land.
!------------------------------------------------------------------
!
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
!
!------------------------------------------------------------------
!     initialize free surface model
!------------------------------------------------------------------
!
      do 300 n=1,3
      do 300 j=1,jmt
      do 300 i=1,imt
        h0(i,j,n) = c0
        u0(i,j,n) = c0
        v0(i,j,n) = c0
300   continue
!
!------------------------------------------------------------------
!     initialize controls
!------------------------------------------------------------------
!
      itt   = 0
      totsec = c0
      totday = c0
      years  = c0
!
      return
      end

      subroutine read_kmt
!
!------------------------------------------------------------------
!     read array 'kmt' from file 'ocean.kmt'.  This is written as:
!     write(56,15)'kmt     ',imt,jmt
!     do j=1,jmt
!     write(56,16)(kmt(i,j),i=1,imt)
!     continue
!------------------------------------------------------------------
!
#include "param.h"
#include "coord.h"
#include "iounit.h"
#include "levind.h"
#include "switch.h"
      character*8 string
!
!  if initial timestep read in kmt field
!
       open(iokmt,file=fnkmt,status='old',iostat = iostat)
       if(iostat.ne.0)then
         write(stdout,11)fnkmt,iostat
         stop
       endif
       read(iokmt,15)string,imt1,jmt1
!
       if(imt1.ne.IMT.or.jmt1.ne.JMT.or.string.ne.'kmt     ')then
         if(string.ne.'kmt     ')then
           write(stdout,12)fnkmt,string
         else
           write(stdout,13)fnkmt,imt1,jmt1,IMT,JMT
         endif
         stop
       endif
       do 20 j=JMT,1,-1
         read(iokmt,16)(kmt(i,j),i=1,IMT)
 20    continue
       close(iokmt)
!
  11  format(1x,'Error opening file "',a,'".',/,'  iostat = ',i4,
     &      '    Program stopping')
  12  format(1x,'Error opening file "',a,'".',/,'  File type faulty',
     &      ' File type = ',a8,/,' Program stopping')
  13  format(1x,'Error opening file "',a,'".',/,'  File dimensions ',
     &      'incompatable with current model.',/,
     &      ' imt, jmt from file =',2i5,/,' imt, jmt of model  =',2i5,/,
     &      ' Program stopping')
  15  format(1x,a8,6i10)
  16  format(1x,30i3)
      end

      subroutine topog
!
!------------------------------------------------------------------
!     construct an "idealized" world ... piece by piece
!     note: this topography will map into arbitrary model
!     resolution
!------------------------------------------------------------------
!
#include "param.h"
#include "coord.h"
#include "grdvar.h"
#include "levind.h"
!
!     define the model topography on the "t" grid using "kmt"
!     set maximum level everywhere except on boundary
!
      do 100 j=2,jmtm1
        do 90 i=2,imtm1
          kmt(i,j) = km
90      continue
100   continue
!
      do 200 j=1,jmt
        kmt(imt,j) = 0
        kmt(1,j) = 0
200   continue
      do 300 i=1,imt
        kmt(i,jmt) = 0
        kmt(i,1) = 0
300   continue
!
!     antarctica
!
      call setkmt (-90.0, 0.0, c360, -80.0, 0.0, c360, 0)
      call setkmt (-80.0, c360-25.0, c360, -70.0, c360, c360, 0)
      call setkmt (-80.0, 0.0, c360, -70.0, 0.0, 170.0, 0)
      call setkmt (-80.0, c360-135.0, c360-60.0, -68.0, c360-75.0,
     &              c360-60.0, 0)
      call setkmt (-70.0, 0.0, 155.0, -67.0, 50.0, 145.0, 0)
!
!     australia
!
      call setkmt (-35.0, 116.0, 120.0, -31.0, 114.0, 130.0, 0)
      call setkmt (-38.0, 140.0, 151.0, -31.0, 130.0, 151.0, 0)
      call setkmt (-31.0, 115.0, 153.0, -20.0, 113.0, 149.0, 0)
      call setkmt (-20.0, 113.0, 149.0, -11.0, 131.0, 143.0, 0)
!
!     south america
!
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
!
!     central america
!
      call setkmt (6.0, c360-78.0, c360-75.0, 20.0, c360-105.0,
     &             c360-97.0, 0)
      call setkmt (20.0, c360-105.0, c360-97.0, 30.0, c360-115.0,
     &             c360-94.0, 0)
!
!     north america
!
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
!
!     greenland
!
      call setkmt (60.0, c360-45.0, c360-45.0, 75.0, c360-58.0,
     &             c360-19.0, 0)
!
!     africa
!
      call setkmt (-35.0, 19.0, 28.0, 6.0, 8.0, 50.0, 0)
      call setkmt (6.0, 0.0, 50.0, 18.0, 0.0, 56.0, 0)
      call setkmt (18.0, 0.0, 56.0, 26.0, 0.0, 59.0, 0)
      call setkmt (6.0, c360-10.0, c360, 18.0, c360-18.0, c360, 0)
      call setkmt (18.0, c360-18.0, c360, 26.0, c360-15.0, c360, 0)
!
!     northern africa &  europe & asia
!
      call setkmt (26.0, c360-15.0, c360, 40.0, c360-7.0, c360, 0)
      call setkmt (40.0, c360-7.0, c360, 50.0, c360, c360, 0)
!
      call setkmt (26.0, 0.0, 126.0, 40.0, 0.0, 122.0, 0)
      call setkmt (40.0, 0.0, 130.0, 50.0, 0.0, 140.0, 0)
      call setkmt (50.0, 0.0, 140.0, 60.0, 8.0, 140.0, 0)
      call setkmt (60.0, 8.0, 163.0, 65.0, 13.0, 180.0, 0)
      call setkmt (65.0, 13.0, 188.0, 70.0, 20.0, 180.0, 0)
      call setkmt (70.0, 70.0, 180.0, 75.0, 90.0, 100.0, 0)
!
!     add an "idealized" ridge in the atlantic & pacific
!
      level = indp (3500.e2, zw, km)
!
      call setkmt (-20.0, c360-20.0, c360-10.0, 30.0, c360-45.0,
     &              c360-35.0, level)
      call setkmt (30.0, c360-45.0, c360-35.0, 60.0, c360-20.0,
     &              c360-30.0, level)
!
      call setkmt (-60.0,c360-100.0, c360-130.0, 40.0, c360-160.0,
     &              180.0, level)
      level = indp (2000.e2, zw, km)
      call setkmt (-50.0, c360-120.0, c360-120.0, 30.0, 190.0, 190.0,
     &             level)
!
!     set cyclic boundary conditions
!
      do 400 j=1,jmt
        kmt(1,j)   = kmt(imtm1,j)
        kmt(imt,j) = kmt(2,j)
400   continue
      return
      end


      function indp(c, ca, ia)
!
!
!==================================================================
!
!     indp = index of nearest data point
!
!     inputs:
!
!     c      = arbitrary data point...same units as ca
!     ca     = array of data points
!              (must be monotonically increasing)
!     ia     = dimension of ca
!
!     output:
!
!     indp =  index of nearest data point to "c"
!             if "c" is outside the domain of "ca" then i
!             ndp is set to 1 or ia depending on  whether
!             c(1) or c(ia) is closest
!
!        note: if "ca" is dimensioned ca(0:ia) in the calling
!              program, then the returned index should be reduced
!              by one to account for the zero base.
!
!     example:
!
!     let model depths be defined by the following:
!     parameter (km=5)
!     dimension z(km)
!     data z /5.0, 10.0, 50.0, 100.0, 250.0/
!
!     k1 = indp (12.5, z, km)
!     k2 = indp (0.0, z, km)
!
!     k1 would be set to 2, & k2 would be set to 1 so that
!     z(k1) would be the nearest data point to 12.5 and z(k2) would
!     be the nearest data point to 0.0
!
!==================================================================
!
!
      dimension ca(ia)
!
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
!
!------------------------------------------------------------------
!     set the topography mask "kmt(i,j)" = "num" within the area of
!     the parallelogram bounded by vertices:
!     (alat1,slon1), (alat1,elon1), (alat2,slon1), & (alat2,elon2)
!------------------------------------------------------------------
!
#include "param.h"
#include "coord.h"
#include "levind.h"
!
!     convert the four vertices into model indices
!     (js,is1), (js,ie1), (je,is2), (je,ie2)
!
      j1 = indp (alat1, yt, jmt)
      j2 = indp (alat2, yt, jmt)
      js = min (j1,j2)
      je = max (j1,j2)
!
      i1  = indp (slon1, xt, imt)
      i2  = indp (elon1, xt, imt)
      is1 = min (i1,i2)
      ie1 = max (i1,i2)
!
      i1  = indp (slon2, xt, imt)
      i2  = indp (elon2, xt, imt)
      is2 = min (i1,i2)
      ie2 = max (i1,i2)
!
      is = is1
      ie = ie1
!
!     fill in the area bounded by (js,is1), (js,ie1),
!                                 (je,is2), (je,ie2)
!
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
!
!==================================================================
!     this subroutine returns estimates of global mean potential
!     temperature for model initialization as a function of depth.
!     it is used to produce a reference thermal stratification for
!     the upper 2000m of the MOM's test case.  below 2000m, the
!     potential temperature returned is 2.0 degrees C.  surface
!     values are set slightly above 18.4 degrees C at the reference
!     latitude "reflat".
!     the estimates are produced from a 7th order ploynomial fit to
!     the annual mean world ocean potential temperature observations
!     of Levitus (1982).
!
!     input [units]:
!       a latitdue (ydeg): [degrees]
!       a zt value (depth): [centimeters]
!     output [units]:
!       potential temperature estimate (est): [degrees centigrade]
!
!     variables:
!       coeft     = coefficients for the polynomial fit of potential
!                   temperature vs. depth
!       reflat    = reference latitude at which observed surface
!                   temperatures approximately equal coeft(1)
!       factor    = the ratio of the cosine of the latitude
!                   requested ("ydeg") to the reference latitude
!                   ("reflat")used to scale the upper 2000 meters
!                   of the vertical temperature profile
!       tmin,tmax = the minumum and maximum potential temperatures
!                   allowed at the time of model initialization
!
!     reference:
!       Levitus, S., Climatological atlas of the world ocean, NOAA
!     Prof. Paper 13, US Gov't printing Office, Washington, DC,
!     1982.
!
      parameter (ndeg=7)
      dimension coeft(ndeg+1)
      save coeft, tmin, tmax, reflat, c0,c1,c2,pi,refcos
      data coeft / 0.184231944E+02,-0.430306621E-01, 0.607121504E-04
     &           ,-0.523806281E-07, 0.272989082E-10,-0.833224666E-14
     &           , 0.136974583E-17,-0.935923382E-22/
      data tmin, tmax, reflat /2.0, 25.0, 34.0/
      data in/1/
!
!==================================================================
!
      if(in.eq.1)then
      in=-1
      c0 = 0.0
      c1 = 1.0
      c2 = 2.0
      pi = atan(1.0) * 4.0
      refcos = abs(cos(pi*reflat/180.))
      endif
!
      coslat = abs(cos(pi*ydeg/180.))
      factor = coslat/refcos
      z = depth * 0.01
!
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
!
      if (est .gt. tmax) est = tmax
      if (est .lt. tmin) est = tmin
!
      theta0 = est
      return
      end


      subroutine tmngr (dt)
!
!
!==================================================================
!
!                      t-i-m-e      m-a-n-a-g-e-r
!
!     inputs:
!
!     itt    = current model time step (from ctmngr.h).
!     dtts   = length of time step in seconds.
!     dgnstc = number of days between printing diagnostics
!     nmix   = number of timesteps between mixint timesteps
!     tsi    = number of days between printing timestep information
!     days   = day at which to end run.
!
!     outputs:
!
!     itt    = itt (on input)  + 1
!     totday = total number of days elapsed
!              (0.5 = noon on first day).
!     years is number of years elapsed (365.5 days = 1 year)
!     diagts = true if time to print diagnostics
!     mixts  = true if a mixing timestep
!     prntsi = true if time to print timestep information
!     last   = true if last timestep of run
!
!==================================================================
!
#include "pconst.h"
#include "ctmngr.h"
#include "switch.h"
#include "timelv.h"
      logical tmnset
!
!------------------------------------------------------------------
!     increment time step counters
!     and create a time stamp
!------------------------------------------------------------------
!
      itt = itt + 1
      totsec = totsec + dt
      totday = totsec*secday
      years  = totday/365.25
      write(stamp,11)totday,years
  11  format('totday =',f8.2,',years =',f8.3)
!
!------------------------------------------------------------------
!
!     set time dependent logical switches.
!     1. set last if it is the last timestep
!     2. set diagts if time to do diagnostics
!     3. set mixts if time mixing to be done
!     4. set prntsi if time to print timestep info
!     5. set snapts if time to archive model
!
!------------------------------------------------------------------
!
      last   = (totday .ge. days)
      diagts = tmnset (dgnstc, totday, dt)
      mixts  = (mod(itt,nmix) .eq. 1)
      prntsi = tmnset (tsi, totday, dt)
      snapts = tmnset (snapd, totday, dt)
      archts = tmnset (archd, totday, dt) .and. .not. last
!
      return
      end


      function tmnset (days, totday, dt)
!
!==================================================================
!     tmnset = true if current time is within 1/2 time step of
!              desired interval
!     days   = desired time interval
!     totday = accumulated time in days from start of run
!     dt     = time step in seconds
!==================================================================
!
#include "pconst.h"
      logical tmnset
!
!  return false if interval days is zero
!
      if(days.eq.c0)then
        tmnset = .false.
        return
      endif
!
!  else return true if the elapsed time 'totday' is within half
!  a timestep 'dt' of an integer multiple of the interval 'days'
!
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



      subroutine snap
!
!==================================================================
!     write snapshot file
!     not currently implemented except for hdf option
!==================================================================
!
#include "iounit.h"
#include "pconst.h"
!
      if(ftsnap.eq.'asc ')then
        write(stdout,*) 'Snapshot write: file type "',
     &                  ftsnap,'" not supported'
#ifdef ncsa_hdf
      elseif(ftsnap.eq.'hdf ')then
        call snap_hdf
#endif
      else
        write(stdout,*) 'Snapshot write: file type "',
     &                  ftsnap,'" not supported'
      endif
      return
      end

      subroutine archrd
!
!==================================================================
!     read archive file
!==================================================================
!
#include "iounit.h"
#include "pconst.h"
!
      if(ftrest.eq.'rst'.or.ftrest.eq.'ufm')then
        call archrd_std
#ifdef ncsa_hdf
      elseif(ftrest.eq.'hdf')then
        call archrd_hdf
#endif
      else
        write(stdout,*) 'Archive read: file type "',
     &                  ftrest,'" not supported'
      endif
      return
      end


      subroutine archwr
!
!==================================================================
!     write archive file
!==================================================================
!
#include "iounit.h"
#include "pconst.h"
!
      if(ftarch.eq.'rst'.or.ftarch.eq.'ufm')then
        call archwr_std
#ifdef ncsa_hdf
      elseif(ftarch.eq.'hdf')then
        call archwr_hdf
#endif
      else
        write(stdout,*) 'Archive write: file type "',
     &                  ftarch,'" not supported'
      endif
      return
      end


      subroutine pdiag
!
!==================================================================
!     write diagnostics
!==================================================================
!
#include "param.h"
#include "coord.h"
#include "ctmngr.h"
#include "frees.h"
#include "timelv.h"
#include "slabs.h"
      dimension ptemp(imt,km)
      character*12 trname(nt)
!
!==================================================================
! based on MOM common "cprnts.h" - only used by this routine
!
      parameter (nlatpr=10)
!
      dimension prlat(nlatpr), prslon(nlatpr), prelon(nlatpr)
     &,         prsdpt(nlatpr), predpt(nlatpr)
      real      slatxy, elatxy, slonxy, elonxy
!
!     default starting & ending coordinates for matrix prints during
!     diagnostic time steps (lats & lons in degrees, depths in cm)
!     see "cprnts.h"
!
      data prlat /-60.0, 0.0, 27.0, 55.0, 6*100.0/
      data prslon /180.0, 180.0, 180.0, 180.0, 6*0.0/
      data prelon /250.0, 250.0, 250.0, 250.0, 6*0.0/
      data prsdpt /nlatpr*0.0/
      data predpt /nlatpr*6000.e2/
      data slatxy, elatxy, slonxy, elonxy /-90.0, 90.0, 3.0, 357.0/
!==================================================================
      data trname/'temperature ','salinity    '/
!
!------------------------------------------------------------------
!       print the solution at specified latitudes
!------------------------------------------------------------------
!
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
!
        do 720 k=1,km
        do 720 i=1,imt
 720    ptemp(i,k) = u(k,i,jj,nc)
        scl = c1
        write (stdout,9100) 'u velocity', itt
     &,   jj, yt(jj), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
        call matrix (ptemp(1,1), imt, is, ie, ks, ke, scl)
!
        do 730 k=1,km
        do 730 i=1,imt
 730    ptemp(i,k) = v(k,i,jj,nc)
        scl = c1
        write (stdout,9100) 'v velocity', itt
     &,   jj, yt(jj), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
        call matrix (ptemp(1,1), imt, is, ie, ks, ke, scl)
      endif
700   continue
!
!-------------------------------------------------------------------
!       print the free surface model solution
!-------------------------------------------------------------------
!
      scl = c100
      is = indp (slonxy, xt, imt)
      ie = indp (elonxy, xt, imt)
      js = indp (slatxy, yt, jmt)
      je = indp (elatxy, yt, jmt)
      write (stdout,8000) ' surface height (m)'
     &, itt, xt(is), xt(ie), yt(js), yt(je), scl
      call matrix (h0(1,1,np0), imt, is, ie, -js, -je, scl)
!
      scl = c1
      write (stdout,8000) ' u velocity (cm/s)'
     &, itt, xt(is), xt(ie), yt(js), yt(je), scl
      call matrix (u0(1,1,np0), imt, is, ie, -js, -je, scl)
!
      write (stdout,8000) ' v velocity (cm/s)'
     &, itt, xt(is), xt(ie), yt(js), yt(je), scl
      call matrix (v0(1,1,np0), imt, is, ie, -js, -je, scl)
!
      return
!
8000  format(1x,a27,1x,'ts=',i7
     &,', lon:',f6.2,' ==> ',f6.2,',   lat:',f6.2,' ==> ',f6.2
     &,', scaling=',1pg10.3)
9000  format(/' Zonal section printouts at specified latitudes:'/)
9100  format(1x,a27,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     &,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     &,', scaling=',1pg10.3)
      end

            subroutine docmnt (nunit, wrdoc, wrc0, wrc1, wrc2)
!
!=======================================================================
!
!     This code arranges data to be written out as a form of
!     documentation summarizing the model characteristics that
!     uniquely define this run.
!     In order to reap the full benefits of this feature, the user will
!     need to review this routine when configuring a model run or adding
!     new options to the code.
!     We envision this feature as being beneficial to those setting up
!     new model runs, comparing different model runs, and for analysis
!     purposes (much of the information used in setting up mapping
!     procedures and other analysis products are contained in the
!     data written by this routine).
!
!     Five namelist writes are available.  Which ones are written
!     is determined from the 4 logical variables passed into the
!     subroutine (wrdoc, wrc0, wrc1, & wrc2).  The unit to which
!     they are written is set by "nunit".
!
!     If "wrdoc" is true, then the "ndocs" & "docums" namelists
!     are written...
!       The 1st, "ndocs" simply contains one variable (ndoc); the
!     dimension of "docum(ndoc)" found in the 2nd namelist write.
!       The 2nd, "docums" consists of a series of 60 character-long
!     lines of information the user should use to list the variables
!     in the subsequent namelist writes, and to describe the model
!     run in general.
!
!     If "wrc0" is true, then the "cntrl0" namelist is written...
!       Namelist "cntrl0" contains the 2 variables needed to restart
!     the model (itt, years).
!
!       The "cntrl1" and "cntrl2" namelists contain information
!     describing the model's configuration and specifications.
!
!     If "wrc1" is true, then the "cntrl1" namelist is written...
!       "cntrl1" contains variables (non-dimensioned) that will be
!     common to most, if not all, incarnations of the model code.
!
!     If "wrc2" is true, then the "cntrl2" namelist is written...
!       "cntrl2" is intended to contain the larger (dimensioned)
!     data that will help identify the model specifics to users
!     and analysis programs alike.
!
!=======================================================================
!  This routine is based on the corresponding mom routine but
!  contains only a few of the standard features
!=======================================================================
!
#include "param.h"
#include "switch.h"
#include "versno.h"
      parameter (ndocum = 16, ifdmax = 50)
!
      character*60 expnam
      character*12 cifdef(ifdmax), blnk12
      logical wrdoc, wrc0, wrc1, wrc2
      data    blnk12 /'            '/
      namelist /conrl0/ itt,years
      namelist /cntrl1/ expnam
!
      expnam = model(1)(1:60)
!
!-----------------------------------------------------------------------
!     list options that are in effect on first timestep
!-----------------------------------------------------------------------
!
      if (first) then
!
        do 200 n=1,ifdmax
          cifdef(n) = blnk12
 200    continue
        nifdef = 0
        cifdef(1) = '0 ifdefs set'
!
#ifdef cyclic
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cyclic      '
#endif
#ifdef simple_world
        nifdef = nifdef + 1
        cifdef(nifdef) = 'simple_world'
#endif
#ifdef oldav
        nifdef = nifdef + 1
        cifdef(nifdef) = 'oldav       '
#endif
#ifdef presetp
        nifdef = nifdef + 1
        cifdef(nifdef) = 'presetp     '
#endif
#ifdef free_eb
        nifdef = nifdef + 1
        cifdef(nifdef) = 'free_eb     '
#endif
#ifdef de_checkbd
        nifdef = nifdef + 1
        cifdef(nifdef) = 'de_checkbd  '
#endif
#ifdef ncsa_hdf
        nifdef = nifdef + 1
        cifdef(nifdef) = 'ncsa_hdf    '
#endif
#ifdef REAL_8
        nifdef = nifdef + 1
        cifdef(nifdef) = 'REAL_8      '
#endif
        write (stdout,'(/,a45)')
     $   ' The following "ifdef" options are in effect:'
        nn = (nifdef+4)/5
        do 300 i=1,nn
          ns = 1+(5*(i-1))
          ne = ns + 4
          if (ne .gt. nifdef) ne = nifdef
          write (stdout,9101) (cifdef(n),n=ns,ne)
300     continue
        write (stdout,'(a21,/)') ' end of options list.'
      endif
      return
9101  format(5(3x,a12))
      end


      subroutine checks
!
!==================================================================
!     do consistency checks
!==================================================================
!
#include "param.h"
#include "scalar.h"
#include "switch.h"
#include "iounit.h"
!
!   check timesteps
!
      if( dtuv.ne.dtts)then
        write (stdout,9000)
     $  '==> Warning: use of unequal time steps implies the transient '
     $, '             response is unimportant and multiple equilibria '
     $, '             do not exist.                                   '
      endif
!
!   check archive and restart file formats
!
      if(ftrest.ne.'rst'.and.ftrest.ne.'ufm'.and.
     &   ftrest.ne.'hdf '.and.init)then
        print *,' ftrest has invalid value'
        print *,' Program stopping... '
        stop
      endif
      if(ftsnap.ne.'hdf'.and.snapd.ne.0.0)then
        print *,' ftsnap has invalid value'
        print *,' Program stopping... '
        stop
      endif
      if(ftarch.ne.'rst'.and.ftarch.ne.'ufm'.and.
     &   ftarch.ne.'hdf '.and.archd.ne.0.0)then
        print *,' ftarch has invalid value'
        print *,' Program stopping... '
        stop
      endif
#ifndef ncsa_hdf
      if(ftrest.eq.'hdf'.and.init)then
        print *,' hdf restart file not supported.  "ncsa_hdf" not set'
        print *,' Program stopping... '
        stop
      endif
      if(ftsnap.eq.'hdf'.and.snapd.ne.0.0)then
        print *,' hdf snapshot file not supported.  "ncsa_hdf" not set'
        print *,' Program stopping... '
        stop
      endif
      if(ftarch.eq.'hdf'.and.(archd.ne.0.0.or.restrt))then
        print *,' hdf restart file not supported.  "ncsa_hdf" not set'
        print *,' Program stopping... '
        stop
      endif
#endif
      return
9000  format (/,(1x,a80))
      end


      subroutine restio
!
!==================================================================
!     read and write restart file
!     NOTE - this version only saves and restores the current
!     timestep.  The model must therefore be restarted with
!     a forward or euler backward timestep.
!==================================================================
!
#include "param.h"
#include "ctmngr.h"
#include "frees.h"
#include "iounit.h"
#include "levind.h"
#include "slabs.h"
#include "switch.h"
#include "timelv.h"
      dimension t1(km,imt,jmt,nt), t2(km,imt,jmt,nt),
     &          t3(km,imt,jmt,nt)
      dimension u1(km,imt,jmt),  u2(km,imt,jmt), u3(km,imt,jmt)
      dimension v1(km,imt,jmt),  v2(km,imt,jmt), v3(km,imt,jmt)
      dimension h01(imt,jmt),    h02(imt,jmt),   h03(imt,jmt)
      dimension u01(imt,jmt),    u02(imt,jmt),   u03(imt,jmt)
      dimension v01(imt,jmt),    v02(imt,jmt),   v03(imt,jmt)
      character*80 file
!
      equivalence (t1(1,1,1,1),t(1,1,1,1,1))
      equivalence (t2(1,1,1,1),t(1,1,1,1,2))
      equivalence (t3(1,1,1,1),t(1,1,1,1,3))
      equivalence (u1(1,1,1),u(1,1,1,1))
      equivalence (u2(1,1,1),u(1,1,1,2))
      equivalence (u3(1,1,1),u(1,1,1,3))
      equivalence (v1(1,1,1),v(1,1,1,1))
      equivalence (v2(1,1,1),v(1,1,1,2))
      equivalence (v3(1,1,1),v(1,1,1,3))
!
      equivalence (h01(1,1),h0(1,1,1))
      equivalence (h02(1,1),h0(1,1,2))
      equivalence (h03(1,1),h0(1,1,3))
      equivalence (u01(1,1),u0(1,1,1))
      equivalence (u02(1,1),u0(1,1,2))
      equivalence (u03(1,1),u0(1,1,3))
      equivalence (v01(1,1),v0(1,1,1))
      equivalence (v02(1,1),v0(1,1,2))
      equivalence (v03(1,1),v0(1,1,3))
!
!------------------------------------------------------------------
!     read restart data set from archive for timesteps itt
!------------------------------------------------------------------
!
      entry archrd_std
!  check np equals 3
      if(np.ne.3)goto 900
!
!------------------------------------------------------------------
!     initialize t,s,u,v.
!          velocities set to zero
!          temperature and salinity set to land values
!          -2.0 and 45.0 ppt
!------------------------------------------------------------------
!
      do 200 n=1,3
      do 200 j=1,jmt
      do 200 i=1,imt
      do 200 k=1,km
        u(k,i,j,n)=c0
        v(k,i,j,n)=c0
        t(k,i,j,1,n)= -2.0
        t(k,i,j,2,n)= 0.01
200   continue
!
!------------------------------------------------------------------
!     initialize free surface model
!------------------------------------------------------------------
!
      do 300 n=1,3
      do 300 j=1,jmt
      do 300 i=1,imt
        h0(i,j,n) = c0
        u0(i,j,n) = c0
        v0(i,j,n) = c0
300   continue
!
!------------------------------------------------------------------
!  read archived data.  Use file fnrest if defined.
!------------------------------------------------------------------
!
      if(fnrest(1:1).ne.' ')then
      open(iorest, file=fnrest, form='UNFORMATTED',
     &       access='SEQUENTIAL')
      else
        open(iorest, form='UNFORMATTED', access='SEQUENTIAL',err=920)
      endif
      if(ftrest.eq.'rst')then
        read(iorest,err=920) itt, totsec, totday, years, kmt
        read(iorest,err=920) h03,  u03,  v03
        read(iorest,err=920) t3,   u3,   v3
      elseif(ftrest.eq.'ufm')then
#ifdef cyclic
        read(iorest,err=920) itt, totsec, totday, years
        read(iorest,err=920) ((kmt(i,j),i=2,imtm1),j=1,jmt)
        read(iorest,err=920) ((h03(i,j),i=2,imtm1),j=1,jmt)
        read(iorest,err=920) ((u03(i,j),i=2,imtm1),j=1,jmt)
        read(iorest,err=920) ((v03(i,j),i=2,imtm1),j=1,jmt)
        do 210 j=1,JMT
 210    read(iorest,err=920)((t(k,i,j,1,3),k=1,km),i=2,imtm1)
        do 220 j=1,JMT
 220    read(iorest,err=920)((t(k,i,j,2,3),k=1,km),i=2,imtm1)
        do 230 j=1,JMT
 230    read(iorest,err=920)((u(k,i,j,3),k=1,km),i=2,imtm1)
        do 240 j=1,JMT
 240    read(iorest,err=920)((v(k,i,j,3),k=1,km),i=2,imtm1)
!
!  set cyclic boundary conditions
!
        do 260 j=1,jmt
          kmt(1,j)   = kmt(imtm1,j)
          kmt(imt,j) = kmt(2,j)
          h03(1,j)   = h03(imtm1,j)
          h03(imt,j) = h03(2,j)
          u03(1,j)   = u03(imtm1,j)
          u03(imt,j) = u03(2,j)
          v03(1,j)   = v03(imtm1,j)
          v03(imt,j) = v03(2,j)
          do 250 k=1,km
            t(k,1,j,1,3)   = t(k,imtm1,j,1,3)
            t(k,imt,j,1,3) = t(k,2,j,1,3)
            t(k,1,j,2,3)   = t(k,imtm1,j,2,3)
            t(k,imt,j,2,3) = t(k,2,j,2,3)
            u(k,1,j,3)   = u(k,imtm1,j,3)
            u(k,imt,j,3) = u(k,2,j,3)
            v(k,1,j,3)   = v(k,imtm1,j,3)
            v(k,imt,j,3) = v(k,2,j,3)
 250      continue
 260    continue
#else
        read(iorest,err=920) itt, totsec, totday, years
        read(iorest,err=920) kmt
        read(iorest,err=920) h03
        read(iorest,err=920) u03
        read(iorest,err=920) v03
        do 210 j=1,JMT
 210    read(iorest,err=920)((t(k,i,j,1,3),k=1,km),i=1,imt)
        do 220 j=1,JMT
 220    read(iorest,err=920)((t(k,i,j,2,3),k=1,km),i=1,imt)
        do 230 j=1,JMT
 230    read(iorest,err=920)((u(k,i,j,3),k=1,km),i=1,imt)
        do 240 j=1,JMT
 240    read(iorest,err=920)((v(k,i,j,3),k=1,km),i=1,imt)
#endif
      endif
      close(iorest)
      write (stdout,21) iorest, itt
      if(mod(itt,nmix).ne.0)goto 900
      return
!
!  error
!
 900  write(stdout,22) np, nc, nm
      stop
 910  write(stdout,23) itt, nmix
      stop
 920  write(stdout,24)
      stop
!
  21  format(1x,/,' Routine restio',/,
     &       ' Restart file read from unit ',i2,'.  timestep =',i7)
  22  format(' Routine restio - error.',/,
     &       ' Trying to restart with np not equal to 3.',/,
     &       ' np, nc, nm  =',3i5,/,' Programme stopping ...')
  23  format(' Routine restio - error',/,
     &       ' Trying to restart without forward timestep,',/,
     &       ' i.e. mod(itt,nmix).ne.0.    itt, nmix =',2i5,/,
     &       ' Programme stopping ...')
  24  format(' Routine restio',/,
     &       ' error while reading restart dataset,',/,
     &       ' Programme stopping ...')
!
!------------------------------------------------------------------
!     write restart data set to archive.
!------------------------------------------------------------------
!
      entry archwr_std
!
!------------------------------------------------------------------
!  create archive file name
!------------------------------------------------------------------
!
      write (file,'(i10.10)')nint(totday)
      file = 'd'//file(7:10)//'.'//ftarch
      open(iorest, file=file, form='UNFORMATTED',
     &       access='SEQUENTIAL',err=950)
!
      write (stdout,31) itt, file, iorest
      if(ftarch.eq.'rst')then
        write(iorest) itt, totsec, totday, years, kmt
        if(np0.eq.1) write(iorest,err=950)  h01, u01, v01
        if(np0.eq.2) write(iorest,err=950)  h02, u02, v02
        if(np0.eq.3) write(iorest,err=950)  h03, u03, v03
        if(np.eq.1)  write(iorest,err=950)  t1, u1, v1
        if(np.eq.2)  write(iorest,err=950)  t2, u2, v2
        if(np.eq.3)  write(iorest,err=950)  t3, u3, v3
      elseif(ftarch.eq.'ufm')then
#ifdef cyclic
        write(iorest) itt, totsec, totday, years
        write(iorest)((kmt(i,j),i=2,imtm1),j=1,jmt)
        if(np0.eq.1) then
          write(iorest,err=950)  ((h01(i,j),i=2,imtm1),j=1,jmt)
          write(iorest,err=950)  ((u01(i,j),i=2,imtm1),j=1,jmt)
          write(iorest,err=950)  ((v01(i,j),i=2,imtm1),j=1,jmt)
        elseif(np0.eq.2)then
          write(iorest,err=950)  ((h02(i,j),i=2,imtm1),j=1,jmt)
          write(iorest,err=950)  ((u02(i,j),i=2,imtm1),j=1,jmt)
          write(iorest,err=950)  ((v02(i,j),i=2,imtm1),j=1,jmt)
        else
          write(iorest,err=950)  ((h03(i,j),i=2,imtm1),j=1,jmt)
          write(iorest,err=950)  ((u03(i,j),i=2,imtm1),j=1,jmt)
          write(iorest,err=950)  ((v03(i,j),i=2,imtm1),j=1,jmt)
        endif
        do 310 j=1,JMT
 310    write(iorest,err=950)((t(k,i,j,1,np),k=1,km),i=2,imtm1)
        do 320 j=1,JMT
 320    write(iorest,err=950)((t(k,i,j,2,np),k=1,km),i=2,imtm1)
        do 330 j=1,JMT
 330    write(iorest,err=950)((u(k,i,j,np),k=1,km),i=2,imtm1)
        do 340 j=1,JMT
 340    write(iorest,err=950)((v(k,i,j,np),k=1,km),i=2,imtm1)
#else
        write(iorest) itt, totsec, totday, years
        write(iorest) kmt
        if(np0.eq.1) then
          write(iorest,err=950)  h01
          write(iorest,err=950)  u01
          write(iorest,err=950)  v01
        elseif(np0.eq.2)then
          write(iorest,err=950)  h02
          write(iorest,err=950)  u02
          write(iorest,err=950)  v02
        else
          write(iorest,err=950)  h03
          write(iorest,err=950)  u03
          write(iorest,err=950)  v03
        endif
        do 310 j=1,JMT
 310    write(iorest,err=950)((t(k,i,j,1,np),k=1,km),i=1,imt)
        do 320 j=1,JMT
 320    write(iorest,err=950)((t(k,i,j,2,np),k=1,km),i=1,imt)
        do 330 j=1,JMT
 330    write(iorest,err=950)((u(k,i,j,np),k=1,km),i=1,imt)
        do 340 j=1,JMT
 340    write(iorest,err=950)((v(k,i,j,np),k=1,km),i=1,imt)
#endif
      endif
      close(iorest, status = 'KEEP')
      return
!
! error
!
 950  write(stdout,32)
      stop
!
  31  format('  Writing restart file for itt = ',i7,/,
     &       '  filename = ',a,'  unit = ',i2)
  32  format(' Routine restio',/,
     &       ' error while writing restart dataset',/,
     &       ' Programme stopping ...')
      end


      subroutine matrix (array, irdim, istrt, im, jstrt, jm, scale)
!
!==================================================================
!
!     matrix is a general two-dimensional array printing routine,
!     where:
!     array = the array to be printed
!     irdim = the 1st dimension of array
!     istrt = the 1st element of the 1st dimension to be printed
!     im    = the last element of the 1st dimension to be printed
!     jstrt = the 1st element of the 2nd dimension to be printed
!     jm    = the last element of the 2nd dimension to be printed
!             the 2nd dimension is printed in reverse order if both
!             jstrt & jm are negative
!     scale = a scaling factor by which array is divided before
!             printing.  (if this is zero, no scaling is done.)
!             if scale=0, 10 columns are printed across in e format
!             if scale>0, 20 columns are printed across in f format
!
!==================================================================
!
#include "param.h"
!
      dimension array(irdim,1000)
!
      if (jstrt*jm .lt. 0) then
        write (stderr,999)  jstrt, jm
        stop '=>matrix'
      endif
!
!     allow for inversion of 2nd dimension
!
      if (jm .lt. 0) then
        js   = -jm
        je   = -jstrt
        jinc = -1
      else
        js   = jstrt
        je   = jm
        jinc = 1
      endif
!
      if (scale .eq. c0) then
        do 100 is=istrt,im,10
          ie = min(is + 9,im)
          write (stdout,9001) (i, i=is,ie)
          do 90 l=js,je,jinc
            write (stdout,9002) l, (array(i,l),i=is,ie)
90        continue
          write (stdout,'(//)')
100     continue
!
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
!
      return
!
999   format (1x,'jstrt=',i5,' jm=',i5,' in matrix')
9001  format(10i13)
9002  format(1x,i2,10(1pe13.5))
9003  format(3x,20i6)
9004  format(1x,i3,1x,20f6.2)
!
      end
