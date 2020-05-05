      program momap
c
c==================================================================
c  MOM's Apple Pie code 
c            - Stripped down(!) version of the MOM (Bryan-Cox-
c              -Semtner) ocean model rewritten for use on an
c              Array Processor.
c            - D J Webb, IOS.
c            - Base Version - 14 May 1992
c @(#) SCCS information - module:  moma.F     
c @(#)                    version: 1.1    date: 05/15/92
c==================================================================
c
#include "param.h"
c
#include "scalar.h"
#include "switch.h"
#include "coord.h"
#include "ctmngr.h"
#include "grdvar.h"
#include "levind.h"
#include "timelv.h"
#include "slabs.h"
c
#include "cdiag.h"
#include "chmix.h"
#include "crelax.h"
#include "cvbc.h"
#include "cvmix.h"
#include "emode.h"
#include "index.h"
#include "iounit.h"
#include "versno.h"
c
      dimension kmz(imt,jmt)
c
      namelist /eddy/   am, ah, fkpm, fkph, cdbot
      namelist /tsteps/ dtts, dtuv, dtsf
      namelist /parms/  mxscan, sor, crit
      namelist /contrl/ init, fnrest, days, restrt, nmix, eb, ncon, 
     &                  tsi,  dgnstc, snaps 
      namelist /iland/ alonis, alatis
c
c     constants usually initialised in block data
c
      model =' Stripped down MOM code for array processor.'//
     &       ' Version 1.0 ... May 1992  '
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
      dtsf = 10800d0
c
      am    = 1.0e9 
      ah    = 2.0e7
      fkph  = 20.0
      fkpm  = 1.0
      cdbot = 0.0
      eb    = .true.
      nmix  = 16
      ncon  = 1
c
c==================================================================
c     begin introductory section
c     open input file and write MOM version information
c==================================================================
c
      open (stdin,file='ocean.in')
      write (stdout,'(/25x,a80/)') model
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
      read  (stdin, parms)
      write (stdout,parms)
      read  (stdin, iland)
      write (stdout,iland)
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
c     timestep, np is next and nm is previous timestep.
c------------------------------------------------------------------
c
      nnm = 1
      nnc = 2
      nnp = 3
      nm = nnm
      nc = nnc
      np = nnp
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
      if (init) then
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
      endif
c
c------------------------------------------------------------------
c     find the island perimeter paths
c     and compute depths and reciprocal depths
c------------------------------------------------------------------
c
      call iperim (kmz)
      do 1400 j=1,jmt
      do 1390 i=1,imt
        hr(i,j) = c0
        h(i,j)  = c0
        if (kmu(i,j) .ne. 0) then
          hr(i,j) = c1/zw(kmu(i,j))
          h (i,j) = zw(kmu(i,j))
        endif
1390  continue
1400  continue
c
c------------------------------------------------------------------
c     compute an array to indicate "interior" stream function grid
c     boxes and set cyclic boundary conditions
c------------------------------------------------------------------
c
      do 1600 i=1,imt
      kmz(i,1)   = 0
      kmz(i,jmt) = 0
1600  continue
      do 1700 j=2,jmtm1
      do 1690 i=2,imu
      kmz(i,j) = min(kmu(i-1,j-1),kmu(i,j-1),kmu(i-1,j),kmu(i,j))
1690  continue
1700  continue
      do 1800 j=1,jmt
      kmz(1,j) = kmz(imtm1,j)
      kmz(imt,j)=kmz(2,j)
1800  continue
c
c-------------------------------------------------------------------
c     compute longitudinal start & end indices for 
c     stream function calculations
c-------------------------------------------------------------------
c
      do 1900 l=1,lseg
        do 1890 j=1,jmt
          isz(j,l) = 0
          iez(j,l) = 0
1890    continue
1900  continue
      lsegp = lseg + 1
c
c     indices for stream function calculation
c
      do 2200 j=3,jscan
        l=1
        do 2190 i=2,imum1
          if (kmz(i-1,j) .eq. 0. .and. kmz(i,j) .ne. 0.) isz(j,l)=i
          if (i .eq. 2 .and. kmz(imt,j) .ne. 0.) isz(j,l)=2
          if (kmz(i,j) .ne. 0. .and. kmz(i+1,j) .eq. 0.) iez(j,l)=i
          if (i .eq. imum1 .and. kmz(i+1,j) .ne. 0.) iez(j,l)=i
          if (kmz(i,j) .ne. 0. .and. kmz(i+1,j) .eq. 0.) l=l+1
          if (l .gt. lsegp) then
            write (stderr,9499) ' stream function',j, lsegp
            stop ' ocean 4'
          endif
2190    continue
2200  continue
c
      write (stdout,9451) ' stream function '
      do 2500 jrev=1,jmt
        j = jmt - jrev + 1
        write (stdout,9461) j, (isz(j,l),iez(j,l),l=1,lseg)
2500  continue
c
c------------------------------------------------------------------
c     initialize various things
c------------------------------------------------------------------
c
      do 2600 j=1,jmt
      do 2590 i=1,imt
      ztd(i,j) = c0
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
        c2dtsf = dtsf
      else
        c2dtts = c2*dtts
        c2dtuv = c2*dtuv
        c2dtsf = c2*dtsf
      endif
c
c------------------------------------------------------------------
c       step model
c------------------------------------------------------------------
c
2800  continue
#ifdef debugprint
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  DEBEG
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print *,'******************************************************'
     &       ,'******************************************************'
      print *,'******************************************************'
     &       ,'******************************************************'
      print * 
      print *,' main routine before call to step:'
      print *,' itt =',itt,' mixts,eb,mxpas2 =',mixts,eb,mxpas2
      print *,' surface b.c. (smf,stf) =',smf,stf
c!!END DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
      print *,' step model.  itt = ',itt
      call step 
#ifdef debugprint
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  DEBEG
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print *
      print *,' main routine following call to step:'
      print *,' itt =',itt,' mixts,eb,mxpas2 =',mixts,eb,mxpas2
      print *,' surface b.c. (smf,stf) =',smf,stf
      print 2801,' np,t(k,60,20,1,np) =',np,(t(k,60,20,1,np),k=1,km)
      print 2801,' nc,t(k,60,20,1,nc) =',nc,(t(k,60,20,1,nc),k=1,km)
      print 2801,' nm,t(k,60,20,1,nm) =',nm,(t(k,60,20,1,nm),k=1,km)
      print 2801,' np,t(k,60,20,2,np) =',np,(t(k,60,20,2,np),k=1,km)
      print 2801,' nc,t(k,60,20,2,nc) =',nc,(t(k,60,20,2,nc),k=1,km)
      print 2801,' nm,t(k,60,20,2,nm) =',nm,(t(k,60,20,2,nm),k=1,km)
      print 2801,' np,p(60,20,np)*e-12=',np,p(60,20,np)*1e-12
      print 2801,' nc,p(60,20,nc)*e-12=',nc,p(60,20,nc)*1e-12
      print 2801,' nm,p(60,20,nm)*e-12=',nm,p(60,20,nm)*1e-12
      print 2801,' np,u(k,60,20,np) =  ',np,(u(k,60,20,np),k=1,km)
      print 2801,' nc,u(k,60,20,nc) =  ',nc,(u(k,60,20,nc),k=1,km)
      print 2801,' nm,u(k,60,20,nm) =  ',nm,(u(k,60,20,nm),k=1,km)
      print 2801,' np,v(k,60,20,np) =  ',np,(v(k,60,20,np),k=1,km)
      print 2801,' nc,v(k,60,20,nc) =  ',nc,(v(k,60,20,nc),k=1,km)
      print 2801,' nm,v(k,60,20,nm) =  ',nm,(v(k,60,20,nm),k=1,km)
2801  format(a,i2,8f12.6,/,23x,8f12.6)
c!!END DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
c
c-----------------------------------------------------------------
c       if this is the end of the 1st pass of an euler backward
c       setup second pass
c------------------------------------------------------------------
c
      if (mixts .and. eb) then
        eots   = .true.
        nc = nnp
        np = nnm
        mixts  = .false.
        mxpas2 = .true.
        go to 2800
      endif
c
c------------------------------------------------------------------
c  at end of second pass of euler backward correct pointers
c------------------------------------------------------------------
c
      if(mxpas2)then
         nc = nnc
         nm = nnp
      endif
c
c-----------------------------------------------------------------
c      print time step information, print diagnostics and archive 
c      as required
c-----------------------------------------------------------------
c
      if (prntsi.and.eots) write (stdout,9602) itt, stamp, ektot,
     &                      dtabs(1), dtabs(2), mscan,
     &                      tddt(1),tddt(2),tvar(1),tvar(2)
      if (diagts.and.eots) call pdiag
      if (snapts.and.eots) call arch
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
c     close all units
c------------------------------------------------------------------
c
      close (unit=stdin,status='keep')
      close (unit=stderr,status='keep')
      write (stdout,'(1x,a)') ' ==> END of model run'
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
9602  format (1x,'ts=',i7, 1x, a32, ', ke=', 1pe13.6,' dtemp=',1pe13.6
     &,        ' dsalt=',1pe13.6,' scans=',i3
     &,/,25x,'tddt(1),tddt(2),tvar(1),tvar(2) =',1p4e13.6)
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
#include "emode.h"
c
#include "cdiag.h"
c  DEBUG
#include "cvbc.h"
c
c-----------------------------------------------------------------
c     1.  calculate barotropic correction for momentum advection 
c         equation
c-----------------------------------------------------------------
c
      do 30 jc=2,jmtm1
      do 30 ic=2,imtm1
      call corv(ic,jc)
30    continue
c
c-----------------------------------------------------------------
c      main loop over grid cells from south to north
c      and from west to east
c      1. set vertical boundary conditions (surface & bottom)
c      2. calculate internal mode velocities
c      3. calculate tracers
c-----------------------------------------------------------------
c
      do 100 jc=2,jmtm1
      do 100 ic=2,imtm1
        call setvbc (ic,jc)
#ifdef debugprint
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  DEBEG
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(ic.eq.20 .and. jc.eq.20)then
      print *
      print *,' subroutine step. after call to setvbc'
      print *,'  ic,jc   =',ic,jc
      print *,'  smf,stf =',smf,stf
      endif
c!!END DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
        call clinic (ic,jc)
        call tracer (ic,jc)
100   continue
c
c-----------------------------------------------------------------
c      calculate stream function 
c-----------------------------------------------------------------
c
      call vort
#ifdef debugprint
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  DEBEG
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print *
      print *,' subroutine step.  following call to vort'
      print 7,' zu(60,20),zv(60,20) =',zu(60,20),zv(60,20)
      print 7,' ztd(60,20)          =',ztd(60,20)
   7  format(a,1p8e16.6)
c!!END DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
      call stream
#ifdef debugprint
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  DEBUG
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print *
      print *,' subroutine step.  following call to stream'
      print *,' nm,nc,np            =',nm,nc,np
      print 7,' p(60,20,nm-nc-np)   =',p(60,20,nm),p(60,20,nc),
     &          p(60,20,np)
c!!END DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
c
c-----------------------------------------------------------------
c     1.  add barotropic velocity to baroclinic 
c-----------------------------------------------------------------
c
      do 300 jc=2,jmtm1
      do 300 ic=2,imtm1
      call addv(ic,jc)
300   continue
c
c-----------------------------------------------------------------
c  set cyclic boundary condition
c-----------------------------------------------------------------
c
      do 420 jc=2,jmt-1
        p(  1,jc,np) = p(imtm1,jc,np)
        p(imt,jc,np) = p(    2,jc,np)
        do 410 k =1,km
          u(k,  1,jc,np) = u(k,imum1,jc,np)
          u(k,imu,jc,np) = u(k,    2,jc,np)
          v(k,  1,jc,np) = v(k,imum1,jc,np)
          v(k,imu,jc,np) = v(k,    2,jc,np)
          do 400 n=1,nt
            t(k,  1,jc,n,np) = t(k,imtm1,jc,n,np)
            t(k,imt,jc,n,np) = t(k,    2,jc,n,np)
400       continue
410     continue
420     continue
c
c-----------------------------------------------------------------
c  collect timestep statistics
c-----------------------------------------------------------------
c
      if(prntsi.and.eots) then
        ektot = c0
        do 500 n=1,nt
        tddt(n)  = c0
        dtabs(n) = c0
        tvar(n)  = c0
 500    continue
        do 510 jc=2,jmtm1
        do 510 ic=2,imtm1
        call diag(ic,jc)
 510    continue
        ektot = ektot/volume
        do 520 n=1,nt
        tddt(n)  = tddt(n)/volume
        dtabs(n) = dtabs(n)/volume
        tvar(n)  = tvar(n)/volume
 520    continue
      endif
c
      return
      end
 
      subroutine clinic(ic,jc)
c
c==================================================================
c
c     clinic computes internal mode velocities at time level n+1
c     for a vertical column of k cells with i and j values equal to
c     ic and jc.
c
c==================================================================
c
c
#include "param.h"
c
#include "scalar.h"
#include "coord.h"
#include "grdvar.h"
#include "levind.h"
#include "timelv.h"
#include "slabs.h"
#include "emode.h"
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
      boxa = dx*csu(jc)*dy
      boxar = c1/boxa
c
c    set up masks, equal to 0 on sub-seabottom points.
c
      do 110 k=1,km
      maskop(k)=0
      if(k.le.kmu(ic,jp))maskop(k)=1
      maskom(k)=0
      if(k.le.kmu(ic,jm))maskom(k)=1
      maskpo(k)=0
      if(k.le.kmu(ip,jc))maskpo(k)=1
      maskmo(k)=0
      if(k.le.kmu(im,jc))maskmo(k)=1
 110  continue 
c
c=============================================================
c  calculate horizontal and vertical advective fluxes
c=============================================================
c
      dxp5 = p5*dx
      dyp5 = p5*dy
      do 160 k=1,kmc
      fvn(k) = (sfv(ic,jp)+(v(k,ic,jc,nc)+v(k,ic,jp,nc))*p5)
     &                               *dx*cst(jp)*dz(k)
      fvs(k) = (sfv(ic,jc)+(v(k,ic,jc,nc)+v(k,ic,jm,nc))*p5)
     &                               *dx*cst(jc)*dz(k)
      fue(k) = (sfu(ip,jc)+(u(k,ic,jc,nc)+u(k,ip,jc,nc))*p5)
     &                               *dy*dz(k)
      fuw(k) = (sfu(ic,jc)+(u(k,ic,jc,nc)+u(k,im,jc,nc))*p5)
     &                               *dy*dz(k)
      fw(k)  = (fue(k)-fuw(k)+fvn(k)-fvs(k))
#ifdef debugprint
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  DEBEG
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(ic.eq.60.and.jc.eq.20.and.k.lt.4)then
      print *
      print *,' subroutine clinic.  after calculating fluxes'
      print *,' k =',k
      print 7,' current at northern boundary'
      print 7,' sfv(ic,jp),v,vop,tot=',sfv(ic,jp),
     &          v(k,ic,jc,nc),v(k,ic,jp,nc),
     &          sfv(ic,jp)+(v(k,ic,jc,nc)+v(k,ic,jp,nc))*p5
      print 7,' current at westhern boundary'
      print 7,' sfu(ic,jc),u,umo,tot=',sfu(ic,jc),
     &          u(k,ic,jc,nc),u(k,im,jc,nc),
     &          sfu(ic,jc)+(u(k,ic,jc,nc)+u(k,im,jc,nc))*p5
      print 7,' fvn(k),fvs(k)       =',fvn(k),fvs(k)
      print 7,' fue(k),fuw(k),fw(k) =',fue(k),fuw(k),fw(k)
   7  format(a,1p8e16.6)
      endif
c!!END DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
 160  continue
c
c=============================================================
c    integrate down from surface and set flux zero
c    at the top (rigid lid) and bottom.
c    note: vertical velocity may be non zero at physical bottom
c=============================================================
c
      fw(0)=c0
      do 180 k=1,km-1
      fw(k)=fw(k)+fw(k-1)
 180  continue
      fw(km)=c0
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
c         compute hydrostatic pressure gradient
c-----------------------------------------------------------------
c
c    1. compute it at the first level
c
      fxa = grav*dzw(0)*csur(jc)
      fxb = grav*dzw(0)*dy2r
      t1  = rhpp(1) - rhoo(1)
      t2  = rhop(1) - rhpo(1)
      dpdx(1) = (t1-t2)*fxa*dx2r
      dpdy(1) = (t1+t2)*fxb
c
c    2. compute the change in pressure gradient between levels
c
      fxa = grav*csur(jc)*p5
      fxb = grav*dy4r
      do 200 k=2,kmc
        t1 = rhpp(k)+rhpp(k-1)- rhoo(k)-rhoo(k-1)
        t2 = rhop(k)+rhop(k-1)- rhpo(k)-rhpo(k-1)
        dpdx(k) = (t1-t2)*fxa*dzw(k-1)*dx2r
        dpdy(k) = (t1+t2)*fxb*dzw(k-1)
 200  continue
c
c=============================================================
c         integrate downward from the first level
c=============================================================
c
      do 220 k=2,kmc
        dpdx(k) = dpdx(k-1) + dpdx(k)
        dpdy(k) = dpdy(k-1) + dpdy(k)
 220  continue
c
c-----------------------------------------------------------------
c         calculate quantities for the computation of
c         vertical diffusion and advection of momentum
c         - positive upwards
c-----------------------------------------------------------------
c
      do 1600 k=1,kmc-1
        vmf(k,1) = fkpm*(u(k+1,ic,jc,nm) - u(k  ,ic,jc,nm))
     &                     *dzwr(k)*dx*csu(jc)*dy
        vmf(k,2) = fkpm*(v(k+1,ic,jc,nm) - v(k  ,ic,jc,nm))
     &                     *dzwr(k)*dx*csu(jc)*dy
        fwb1(k) = fw(k)*(u(k,ic,jc,nc) + u(k+1,ic,jc,nc))*p5
        fwb2(k) = fw(k)*(v(k,ic,jc,nc) + v(k+1,ic,jc,nc))*p5
1600  continue
c
c     set k=0 elements to reflect wind stress
c     set level "kz" for a no bottom drag condition.
c     set k=0 elements of vertical advection array to the rigid lid
c     condition (w(0) = 0). also assume w(km)=0. If "kz" < "km"
c     and there is a bottom slope, w(kz) may be non zero.
c     smf is the surface stress acting on the ocean
c     bmf is the bottom stress acting on the ocean floor
c
      vmf(0,1)  = -smf(1)*boxa
      vmf(0,2)  = -smf(2)*boxa
      vmf(kmc,1)= -bmf(1)*boxa
      vmf(kmc,2)= -bmf(2)*boxa
      fwb1(0)  = c0
      fwb2(0)  = c0
      fwb1(km) = c0
      fwb2(km) = c0
c
c=============================================================
c    main calculation to timestep baroclinic velocity
c=============================================================
c
      do 400 k=1,kmc
      u(k,ic,jc,np) = 
c  coriolis term and pressure gradient
     &   fcor(jc)*v(k,ic,jc,nc) - dpdx(k)
c  u-advection
     & +(fuw(k)*(u(k,im,jc,nc) + u(k,ic,jc,nc))*maskmo(k)*p5
     & - fue(k)*(u(k,ip,jc,nc) + u(k,ic,jc,nc))*maskpo(k)*p5
c  v-advection
     & + fvs(k)*(u(k,ic,jm,nc) + u(k,ic,jc,nc))*maskom(k)*p5
     & - fvn(k)*(u(k,ic,jp,nc) + u(k,ic,jc,nc))*maskop(k)*p5
c  w-advection
     & + fwb1(k  ) - fwb1(k-1)
c  x diffusion  
     & + bbu(jc)*(u(k,ip,jc,nm)-u(k,ic,jc,nm))*maskpo(k)
     & + bbu(jc)*(u(k,im,jc,nm)-u(k,ic,jc,nm))*maskmo(k)
c  y diffusion
     & + ccu(jc)*(u(k,ic,jp,nm)-u(k,ic,jc,nm))*maskop(k)
     & + ddu(jc)*(u(k,ic,jm,nm)-u(k,ic,jc,nm))*maskom(k)
c  z diffusion  
     &  + vmf(k  ,1) - vmf(k-1,1) )*boxar*dzr(k)
c
c***************************************************************
c  the terms that allow for the curvature of the grid near the
c  poles have been left out of this version of the model
c***************************************************************
c
 400  continue
c
      do 500 k=1,kmc
      v(k,ic,jc,np) = 
c  coriolis term and pressure gradient
     &  -fcor(jc)*u(k,ic,jc,nc) - dpdy(k)
c  u-advection
     & +(fuw(k)*(v(k,im,jc,nc) + v(k,ic,jc,nc))*maskmo(k)*p5
     & - fue(k)*(v(k,ip,jc,nc) + v(k,ic,jc,nc))*maskpo(k)*p5
c  v-advection
     & + fvs(k)*(v(k,ic,jm,nc) + v(k,ic,jc,nc))*maskom(k)*p5
     & - fvn(k)*(v(k,ic,jp,nc) + v(k,ic,jc,nc))*maskop(k)*p5
c  w-advection
     & + fwb2(k  ) - fwb2(k-1)
c  x diffusion  
     & +(bbu(jc)*(v(k,ip,jc,nm)-v(k,ic,jc,nm))*maskpo(k)
     & + bbu(jc)*(v(k,im,jc,nm)-v(k,ic,jc,nm))*maskmo(k)
c  y diffusion
     & + ccu(jc)*(v(k,ic,jp,nm)-v(k,ic,jc,nm))*maskop(k)
     & + ddu(jc)*(v(k,ic,jm,nm)-v(k,ic,jc,nm))*maskom(k))
     & *dz(k)
c  z diffusion
     & + vmf(k  ,2) - vmf(k-1,2) )*boxar*dzr(k)
c
c***************************************************************
c  the terms that allow for the curvature of the grid near the
c  poles have been left out of this version of the model
c***************************************************************
c
#ifdef debugprint
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c DEBUG
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(ic.eq.60.and.jc.eq.20.and.k.lt.4)then
      boxvol=boxa*dz(k)
      boxvr=boxar*dzr(k)
      print *
      print *,' subroutine clinic.  main calculation'
      print *,'  ic,jc,k=',ic,jc,k
      print 7,' dx,dy,cst(jc),boxa,boxvol',
     &         dx,dy,cst(jc),boxa,boxvol,
     &         boxvr,boxar*dzr(k)
      print 7,' rhoo,rhpo,rhop,rhpp;k; =',rhoo(k),
     &         rhpo(k),rhop(k),rhpp(k)
      print 7,' fcor*boxvol,dpdy*boxvol (+/boxvol)',
     &- fcor(jc)*u(k,ic,jc,nc)*boxvol, - dpdy(k)*boxvol,
     &- fcor(jc)*u(k,ic,jc,nc),        - dpdy(k)
c  u-advection
      print 7,' u-advc',
     & + fuw(k)*(v(k,im,jc,nc) + v(k,ic,jc,nc))*maskmo(k)*p5,
     & - fue(k)*(v(k,ip,jc,nc) + v(k,ic,jc,nc))*maskpo(k)*p5,
     & +fuw(k)*(v(k,im,jc,nc)+v(k,ic,jc,nc))*maskmo(k)*p5*boxvr,
     & -fue(k)*(v(k,ip,jc,nc)+v(k,ic,jc,nc))*maskpo(k)*p5*boxvr
c  v-advection
      print 7,' v-advc',
     & + fvs(k)*(v(k,ic,jm,nc) + v(k,ic,jc,nc))*maskom(k)*p5,
     & - fvn(k)*(v(k,ic,jp,nc) + v(k,ic,jc,nc))*maskop(k)*p5,
     & +fvs(k)*(v(k,ic,jm,nc)+v(k,ic,jc,nc))*maskom(k)*p5*boxvr,
     & -fvn(k)*(v(k,ic,jp,nc)+v(k,ic,jc,nc))*maskop(k)*p5*boxvr
c  w-advection
      print 7,' w-advc',
     & +fwb2(k),-fwb2(k-1),+fwb2(k)*boxvr,-fwb2(k-1)*boxvr
c  x diffusion  
      print 7,' x-diff',
     & + bbu(jc)*(v(k,ip,jc,nm)-v(k,ic,jc,nm))*maskpo(k),
     & + bbu(jc)*(v(k,im,jc,nm)-v(k,ic,jc,nm))*maskmo(k),
     & bbu(jc)*(v(k,ip,jc,nm)-v(k,ic,jc,nm))*maskpo(k),
     & bbu(jc)*(v(k,im,jc,nm)-v(k,ic,jc,nm))*maskmo(k)*boxvr
c  y diffusion
      print 7,' y-diff',
     & + ccu(jc)*(v(k,ic,jp,nm)-v(k,ic,jc,nm))*maskop(k),
     & + ddu(jc)*(v(k,ic,jm,nm)-v(k,ic,jc,nm))*maskom(k),
     & ccu(jc)*(v(k,ic,jp,nm)-v(k,ic,jc,nm))*maskop(k)*boxvr,
     & ddu(jc)*(v(k,ic,jm,nm)-v(k,ic,jc,nm))*maskom(k)*boxvr
c  z diffusion
      print 7,' z-diff',
     & + vmf(k,2), - vmf(k-1,2),
     & + vmf(k,2)*boxvr, - vmf(k-1,2)*boxvr
      endif
c!!END DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
 500  continue
c
c-----------------------------------------------------------------
c     form time change of vertically averaged forcing
c-----------------------------------------------------------------
c
      zu(ic,jc) = c0
      zv(ic,jc) = c0
      do 600 k=1,kmc
      fx = c2dtsf*dz(k)
      zu(ic,jc) = zu(ic,jc) + u(k,ic,jc,np)*fx
      zv(ic,jc) = zv(ic,jc) + v(k,ic,jc,np)*fx
 600  continue
      zu(ic,jc) = zu(ic,jc)*hr(ic,jc)
      zv(ic,jc) = zv(ic,jc)*hr(ic,jc)
c
c-----------------------------------------------------------------
c     timestep baroclinit terms
c-----------------------------------------------------------------
c
      do 700 k=1,kmc
      q1temp = u(k,ic,jc,np)
      q2temp = v(k,ic,jc,np)
      u(k,ic,jc,np) = u(k,ic,jc,nm) + c2dtuv*u(k,ic,jc,np) 
      v(k,ic,jc,np) = v(k,ic,jc,nm) + c2dtuv*v(k,ic,jc,np) 
#ifdef debugprint
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c DEBUG
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(ic.eq.60.and.jc.eq.20.and.k.lt.4)then
      print *
      print *,' subroutine clinic. before removing vertical mean'
      print *,'  ic,jc,k=',ic,jc,k
      print 7,'  u(k,ic,jc,nm), v(k,ic,jc,nm)     =',
     &           u(k,ic,jc,nm), v(k,ic,jc,nm)
      print 7,'  u(k,ic,jc,np), v(k,ic,jc,np) old =',q1temp,q2temp
      print 7,'  c2dtuv                           =',c2dtuv
      print 7,'  u(k,ic,jc,np), v(k,ic,jc,np) new =',
     &           u(k,ic,jc,np), v(k,ic,jc,np)

      print 7,' zu(ic,jc), zv(ic,jc)              =',
     &          zu(ic,jc),zv(ic,jc)
      endif
c!!END DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
 700  continue
c
c-----------------------------------------------------------------
c     remove vertical mean
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
      q1temp = u(k,ic,jc,np)
      q2temp = v(k,ic,jc,np)
      u(k,ic,jc,np) = u(k,ic,jc,np) - sumu
      v(k,ic,jc,np) = v(k,ic,jc,np) - sumv
#ifdef debugprint
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c DEBUG
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(ic.eq.60.and.jc.eq.20.and.k.lt.4)then
      print *
      print *,' subroutine clinic.  after removing vertical mean'
      print *,'  ic,jc,k=',ic,jc,k
      print 7,'  u(k,ic,jc,np), v(k,ic,jc,np) old =',
     &           q1temp,q2temp
      print 7,'  sumu,sumv                        =',sumu,sumv
      print 7,'  u(k,ic,jc,np), v(k,ic,jc,np) new =',
     &           u(k,ic,jc,np), v(k,ic,jc,np)
      endif
c!!END DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
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
#ifdef debugprint
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  DEBUG
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      data in/1/
      if(ic.eq.60.and.jc.eq.20)then
      print *
      print *,' subroutine tracer, on entry'
      print *,'  ic,jc =',ic,jc
      endif
c!!END DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
c
c    outer loops over i and j indices
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
c    set up masks, equal to 0 on sub-seabottom points.
c
      do 110 k=1,km
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
      fw(k)  = (fue(k)-fuw(k)+fvn(k)-fvs(k))
 160  continue
c
c=============================================================
c    integrate down from surface and set flux zero
c    at the top (rigid lid) and bottom.
c=============================================================
c
      fw(0)=c0
      do 180 k=1,kmc-1
      fw(k)=fw(k)+fw(k-1)
 180  continue
      fw(kmc)=c0
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
        fwb1(k) = fw(k)*(t(k,ic,jc,m,nm)+t(k+1,ic,jc,m,nm))*p5
 930  continue
c
c=============================================================
c  surface and bottom advective and diffusive boundary conditions
c  stf is flux into the ocean
c=============================================================
c
      vtf(0)   = -stf(m)*boxa
      vtf(kmc) = c0
      fwb1(0)  = c0
      fwb1(kmc)= c0
c
c=============================================================
c  main timestepping calculation
c=============================================================
c
      do 960 k=1,kmc
      t(k,ic,jc,m,np) = t(k,ic,jc,m,nm)
     &                + c2dtts*boxar*dzr(k)*(
c  u-advection
     &   fuw(k)*(t(k,im,jc,m,nc) + t(k,ic,jc,m,nc))*maskmo(k)*p5
     & - fue(k)*(t(k,ip,jc,m,nc) + t(k,ic,jc,m,nc))*maskpo(k)*p5
c  v-advection
     & + fvs(k)*(t(k,ic,jm,m,nc) + t(k,ic,jc,m,nc))*maskom(k)*p5
     & - fvn(k)*(t(k,ic,jp,m,nc) + t(k,ic,jc,m,nc))*maskop(k)*p5
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
#ifdef debugprint
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c DEBUG
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if((ic.eq.60.and.jc.eq.20.and.k.lt.4)
c     & .or.   (ic.eq.73.and.jc.eq.5)
c     & .or.   (m.eq.2.and.(t(k,ic,jc,m,np).ne.t(k,ic,jc,m,nm)))
     &      )then
      in=in+1
      boxvr=boxar*dzr(k)
      print *
      print *,' subroutine tracer, main calculation'
      print *,'  ic,jc,k,m =',ic,jc,k,m
      print 6,'  tp, tm    =',t(k,ic,jc,m,np), t(k,ic,jc,m,nm), c2dtts,
     &        boxa*dz(k),boxar*dzr(k)
c  u-advection
      print 7,'  u-advc',
     &   fuw(k)*(t(k,im,jc,m,nc) + t(k,ic,jc,m,nc))*maskmo(k)*p5
     &,- fue(k)*(t(k,ip,jc,m,nc) + t(k,ic,jc,m,nc))*maskpo(k)*p5
     &,  fuw(k)*(t(k,im,jc,m,nc) + t(k,ic,jc,m,nc))*maskmo(k)*p5
     & - fue(k)*(t(k,ip,jc,m,nc) + t(k,ic,jc,m,nc))*maskpo(k)*p5
      print 7,'  /vol  ',
     &   fuw(k)*(t(k,im,jc,m,nc) + t(k,ic,jc,m,nc))*maskmo(k)*p5
     &   *boxvr
     &,- fue(k)*(t(k,ip,jc,m,nc) + t(k,ic,jc,m,nc))*maskpo(k)*p5
     &   *boxvr
     &,( fuw(k)*(t(k,im,jc,m,nc) + t(k,ic,jc,m,nc))*maskmo(k)*p5
     & - fue(k)*(t(k,ip,jc,m,nc) + t(k,ic,jc,m,nc))*maskpo(k)*p5)
     &   *boxvr
   7  format(a,1p8e18.8)
c  v-advection
      print 7,'  v-advc',
     & + fvs(k)*(t(k,ic,jm,m,nc) + t(k,ic,jc,m,nc))*maskom(k)*p5
     &,- fvn(k)*(t(k,ic,jp,m,nc) + t(k,ic,jc,m,nc))*maskop(k)*p5
     &,+ fvs(k)*(t(k,ic,jm,m,nc) + t(k,ic,jc,m,nc))*maskom(k)*p5
     & - fvn(k)*(t(k,ic,jp,m,nc) + t(k,ic,jc,m,nc))*maskop(k)*p5
      print 7,'  /vol  ',
     & + fvs(k)*(t(k,ic,jm,m,nc) + t(k,ic,jc,m,nc))*maskom(k)*p5
     &   *boxvr
     &,- fvn(k)*(t(k,ic,jp,m,nc) + t(k,ic,jc,m,nc))*maskop(k)*p5
     &   *boxvr
     &,( fvs(k)*(t(k,ic,jm,m,nc) + t(k,ic,jc,m,nc))*maskom(k)*p5
     & - fvn(k)*(t(k,ic,jp,m,nc) + t(k,ic,jc,m,nc))*maskop(k)*p5)
     &   *boxvr
c  w-advection
      print 7,'  w-advc',
     & + fwb1(k  ), - fwb1(k-1)
     &,+ fwb1(k  )  - fwb1(k-1)
      print 7,'  /vol  ',
     & + fwb1(k  )*boxvr, - fwb1(k-1)*boxvr
     &,+ (fwb1(k  )  - fwb1(k-1))*boxvr
c  x diffusion
      print 7,'  x-diff',
     & + bbt(jc)*(t(k,ip,jc,m,nm)-t(k,ic,jc,m,nm))*maskpo(k)*dz(k),
     & + bbt(jc)*(t(k,im,jc,m,nm)-t(k,ic,jc,m,nm))*maskmo(k)*dz(k),
     & + bbt(jc)*(t(k,ip,jc,m,nm)-t(k,ic,jc,m,nm))*maskpo(k)*dz(k)
     & + bbt(jc)*(t(k,im,jc,m,nm)-t(k,ic,jc,m,nm))*maskmo(k)*dz(k)
      print 7,'  /vol  ',
     & + bbt(jc)*(t(k,ip,jc,m,nm)-t(k,ic,jc,m,nm))*maskpo(k)*dz(k)
     &   *boxvr
     &,  bbt(jc)*(t(k,im,jc,m,nm)-t(k,ic,jc,m,nm))*maskmo(k)*dz(k)
     &   *boxvr
     &,+(bbt(jc)*(t(k,ip,jc,m,nm)-t(k,ic,jc,m,nm))*maskpo(k)
     & + bbt(jc)*(t(k,im,jc,m,nm)-t(k,ic,jc,m,nm))*maskmo(k))
     &   *boxvr*dz(k)
c  y diffusion
      print 7,'  y-diff',
     & + cct(jc)*(t(k,ic,jp,m,nm)-t(k,ic,jc,m,nm))*maskop(k)*dz(k),
     & + ddt(jc)*(t(k,ic,jm,m,nm)-t(k,ic,jc,m,nm))*maskom(k)*dz(k),
     & + cct(jc)*(t(k,ic,jp,m,nm)-t(k,ic,jc,m,nm))*maskop(k)*dz(k)
     & + ddt(jc)*(t(k,ic,jm,m,nm)-t(k,ic,jc,m,nm))*maskom(k)*dz(k)
      print 7,'  /vol  ',
     & + cct(jc)*(t(k,ic,jp,m,nm)-t(k,ic,jc,m,nm))*maskop(k)
     &   *boxvr*dz(k)
     &,  ddt(jc)*(t(k,ic,jm,m,nm)-t(k,ic,jc,m,nm))*maskom(k)
     &   *boxvr*dz(k)
     &, (cct(jc)*(t(k,ic,jp,m,nm)-t(k,ic,jc,m,nm))*maskop(k)
     & + ddt(jc)*(t(k,ic,jm,m,nm)-t(k,ic,jc,m,nm))*maskom(k))
     &   *boxvr*dz(k)
c  z diffusion  
      print 7,'  z-diff',
     & + vtf(k  ), - vtf(k-1),
     & + vtf(k  )  - vtf(k-1)
      print 7,'  /vol  ',
     & + vtf(k  )*boxvr, - vtf(k-1)*boxvr,
     & + (vtf(k  )  - vtf(k-1))*boxvr
c  composite terms
      temp1=
c  u-advection
     &   fuw(k)*(t(k,im,jc,m,nc) + t(k,ic,jc,m,nc))*maskmo(k)*p5
     & - fue(k)*(t(k,ip,jc,m,nc) + t(k,ic,jc,m,nc))*maskpo(k)*p5
c  v-advection
     & + fvs(k)*(t(k,ic,jm,m,nc) + t(k,ic,jc,m,nc))*maskom(k)*p5
     & - fvn(k)*(t(k,ic,jp,m,nc) + t(k,ic,jc,m,nc))*maskop(k)*p5
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
      temp2 = max(
c  u-advection
     & abs(fuw(k)*(t(k,im,jc,m,nc) + t(k,ic,jc,m,nc))*maskmo(k)*p5),
     & abs(fue(k)*(t(k,ip,jc,m,nc) + t(k,ic,jc,m,nc))*maskpo(k)*p5),
c  v-advection
     & abs(fvs(k)*(t(k,ic,jm,m,nc) + t(k,ic,jc,m,nc))*maskom(k)*p5),
     & abs(fvn(k)*(t(k,ic,jp,m,nc) + t(k,ic,jc,m,nc))*maskop(k)*p5),
c  w-advection
     & abs(fwb1(k  )), abs(fwb1(k-1)),abs(
c  x diffusion
     & +(bbt(jc)*(t(k,ip,jc,m,nm)-t(k,ic,jc,m,nm))*maskpo(k)
     & + bbt(jc)*(t(k,im,jc,m,nm)-t(k,ic,jc,m,nm))*maskmo(k)
c  y diffusion
     & + cct(jc)*(t(k,ic,jp,m,nm)-t(k,ic,jc,m,nm))*maskop(k)
     & + ddt(jc)*(t(k,ic,jm,m,nm)-t(k,ic,jc,m,nm))*maskom(k))
     & * dz(k)
c  z diffusion  
     & + vtf(k  ) - vtf(k-1))   )
      print *,'  nm,nc,np = ',nm,nc,np
      print *, 'tp-tm, temp1*c2dtts*boxvr,temp2*c2dtts*boxvr'
     &,t(k,ic,jc,m,np) - t(k,ic,jc,m,nm)
     &,temp1*c2dtts*boxvr,temp2*c2dtts*boxvr
      print *, ' temp1, temp2'
     &,          temp1, temp2
        if(in.gt.100)then
          print *,' in =',in,' stopping'
          stop
        endif
      endif
c!!END DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
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
#ifdef debugprint
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c DEBUG
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c DEBUG !!!!!! set m to 1 or 2
      if(m.eq.0)then
        in=in+1
        print *,' subroutine tracer - convection'
        print *,' n,ks,m,k,ic,jc = ',n,ks,m,k,ic,jc
        print 6,' den(l),den(k+1), diff =',temp(k),temp(k+1),
     &            temp(k)-temp(k+1)
        print 6,' t(k),  t(k+1)  , diff =',t(k,ic,jc,1,np),
     &     t(k+1,ic,jc,1,np),t(k,ic,jc,1,np)-t(k+1,ic,jc,1,np)
        print 6,' s(k),  s(k+1)  , diff =',t(k,ic,jc,2,np),
     &     t(k+1,ic,jc,2,np),t(k,ic,jc,2,np)-t(k+1,ic,jc,2,np)
   6  format(a,1p3e25.17,(/,24x,3e25.17))
        if(in.gt.200)stop
      endif
c!!END DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
          t(k,ic,jc,m,np) = (dz(k)*t(k,ic,jc,m,np)+
     &                    +  dz(k+1)*t(k+1,ic,jc,m,np))*dzw2r(k)
          t(k+1,ic,jc,m,np) = t(k,ic,jc,m,np)
        endif
 1070   continue
      endif
c
c=============================================================
c  finished !
c=============================================================
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
c     so... this routine should be used only for horizontal 
c     gradients of density.
c
c     inputs:
c
c     t = the input column of temperatures
c     s = the input column of salinities (units: (ppt-35)/1000)
c     ind = 0 for normal density calculations
c         = 1 for comparing levels 1 to 2, 3 to 4, etc.
c         = 2 for comparing levels 2 to 3, 4 to 5, etc.
c           (if ind = 1 or 2, the coefficients for the lower
c            of the 2 levels are used)
c     kk = bottom level to be calculated
c
c     note: salinities are converted to (s-35.0)/1000 for
c     use with constants given in "dncoef.h".  The latter are
c     assumed to use the normal GFDL model salinity units.
c
c     output:
c
c     rho = density with depth dependent offset
c
c=============================================================
c
c
#include "param.h"
      dimension t(km),  s(km),  rho(km)
      dimension to(km), so(km), c(km,9),
     &          to1(km),so1(km),co1(km,9),
     &          to2(km),so2(km),co2(km,9)
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
        sq=(s(k)-35.0)/1000.0-so(k)
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
        sq=(s(k)-35.0)/1000.0-so1(k)
        rho(k) = (co1(k,1) + (co1(k,4)    + co1(k,7)*sq)*sq
     &         + (co1(k,3) +  co1(k,8)*sq + co1(k,6)*tq)*tq)*tq
     &         + (co1(k,2) + (co1(k,5)    + co1(k,9)*sq)*sq)*sq
 500    continue
c
c  compare 2 with 3, 4 with 5 etc
      elseif(ind.eq.2.and.kk.gt.2)then
        do 600 k1=1,kk-2,2
        do 600 k2=1,2
        k=k1+k2
        tq=t(k)-to2(k)
        sq=(s(k)-35.0)/1000.0-so2(k)
        rho(k) = (co2(k,1) + (co2(k,4)    + co2(k,7)*sq)*sq
     &         + (co2(k,3) +  co2(k,8)*sq + co2(k,6)*tq)*tq)*tq
     &         + (co2(k,2) + (co2(k,5)    + co2(k,9)*sq)*sq)*sq
 600    continue
      endif
      return
   99 format(/,' Routine state. Variables "ind" or "kk" out',
     &       ' of range',/,' ind =',i10,' kk =',i10)
      end

      subroutine vort
c
c==================================================================
c
c     compute vorticity for relaxation
c
c==================================================================
c
c
#include "param.h"
c
#include "grdvar.h"
#include "emode.h"
c
c------------------------------------------------------------------
c       ensure zu and zv are cyclic
c------------------------------------------------------------------
c
      do 50 j=2,jmt-1
      zu(1  ,j) = zu(imtm1,j)
      zu(imt,j) = zu(    2,j)
      zv(1  ,j) = zv(imtm1,j)
      zv(imt,j) = zv(   m2,j)
  50  continue
c
c------------------------------------------------------------------
c       form curl of time change in vertically averaged equations
c------------------------------------------------------------------
c
c       all vorticity pts. are computed so that those needed for 
c       the line integral of hole relaxation (immediately adjacent
c       to will be defined).
c
      do 200 j=2,jmtm1
      do 100 i=2,imtm1
      ztd(i,j) = (zu(i,j  ) + zu(i-1,j  ))*dx*csu(j  )
     &          -(zu(i,j-1) + zu(i-1,j-1))*dx*csu(j-1)
      ztd(i,j) =((zv(i,j  ) - zv(i-1,j  )
     &           +zv(i,j-1) - zv(i-1,j-1))*dy
     &           -ztd(i,j))*dx2r*dyr*cstr(j)
100   continue
200   continue
      return
      end
      subroutine stream
c
c==================================================================
c
c     stream is similar to relax except it solves the equation 
c     along diagonals to allow for complete vectorization & 
c     parallelization
c
c==================================================================
c
c
#include "param.h"
c
#include "ctmngr.h"
#include "scalar.h"
#include "switch.h"
#include "grdvar.h"
#include "timelv.h"
c
#include "crelax.h"
#include "emode.h"
#include "index.h"
c
c==================================================================
c     begin introductory section to prepare for the relaxation
c==================================================================
c
      if (first .or. mixts .or. mod(itt,nmix) .eq. 2) then
c
c------------------------------------------------------------------
c       initialize the work area
c------------------------------------------------------------------
c
        do 100 j=1,jmt
          do 90 i=1,imt
            cfn(i,j)  = c0
            cfs(i,j)  = c0
            cfe(i,j)  = c0
            cfw(i,j)  = c0
            rmsk(i,j) = c0
            cof(i,j)  = c1
90        continue
100     continue
c
c------------------------------------------------------------------
c       generate arrays of coefficients for relaxation
c------------------------------------------------------------------
c
c       1st, compute coefficients of the laplacian star
c       (hold non-interior points to zero using start and end 
c       indices)
c
        do 600 j=3,jmtm1
          do 590 l=1,lseg
            if (isz(j,l) .ne. 0) then
              ie = iez(j,l)
              fxa = c2*cstr(j)*cstr(j)
              fxb = c2*csu(j  )*cstr(j)*dyr*dyr
              fxc = c2*csu(j-1)*cstr(j)*dyr*dyr
              do 560 i=isz(j,l),iez(j,l)
                rmsk(i,j) = c1
                cfn(i,j)  = fxb/(h(i-1,j)+h(i,j))
                cfs(i,j)  = fxc/(h(i-1,j-1)+h(i,j-1))
                cfe(i,j)  = fxa*dxr*dxr/(h(i,j)+h(i,j-1))
                cfw(i,j)  = fxa*dxr*dxr/(h(i-1,j)+h(i-1,j-1))
                cpf(i,j) = c1/(cfn(i,j)+cfs(i,j)+cfe(i,j)+cfw(i,j))
560          continue
c
c             normalize coefficients
c
              do 580 i=isz(j,l),iez(j,l)
                cfn(i,j) = cfn(i,j)*cpf(i,j)
                cfs(i,j) = cfs(i,j)*cpf(i,j)
                cfe(i,j) = cfe(i,j)*cpf(i,j)
                cfw(i,j) = cfw(i,j)*cpf(i,j)
580           continue
            endif
590       continue
600     continue
c
c       modify coefficients to include island perimeter points
c
        do 700 isle=1,nisle
          cofis(isle) = c0
          do 690 n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            fxa = c2*cstr(j)*cstr(j)
            fxb = c2*csu(j  )*dyr*dyr*cstr(j)
            fxc = c2*csu(j-1)*dyr*dyr*cstr(j)
            fxd = c0
c           fxd = -c2dtsf*acor*cstr(j)*dyr*omega
            if(hr(i-1,j  ) .ne. c0 .or. hr(i  ,j  ) .ne. c0)
     &        cfn(i,j) = fxb/(h(i-1,j)+h(i,j))
            if(hr(i-1,j-1) .ne. c0 .or. hr(i  ,j-1) .ne. c0)
     &        cfs(i,j) = fxc/(h(i-1,j-1)+h(i,j-1))
            if(hr(i  ,j  ) .ne. c0 .or. hr(i  ,j-1) .ne. c0)
     &        cfe(i,j) = fxa*dxr*dxr/(h(i,j)+h(i,j-1))
            if(hr(i-1,j  ) .ne. c0 .or. hr(i-1,j-1) .ne. c0)
     &        cfw(i,j) = fxa*dxr*dxr/(h(i-1,j)+h(i-1,j-1))
c
            cpf(i,j) = c1/(cfn(i,j)+cfs(i,j)+cfe(i,j)+cfw(i,j))
c
            cfn(i,j) = cfn(i,j)*cpf(i,j)
            cfs(i,j) = cfs(i,j)*cpf(i,j)
            cfe(i,j) = cfe(i,j)*cpf(i,j)
            cfw(i,j) = cfw(i,j)*cpf(i,j)
            cof(i,j) = cst(j)*dx*dy/cpf(i,j)
            cofis(isle) = cofis(isle)+cof(i,j)
690       continue
          cofis(isle) = c1/cofis(isle)
700     continue
      endif
c
c------------------------------------------------------------------
c     normalize vorticity "ztd"
c------------------------------------------------------------------
c
      do 800 j=1,jmt
        do 790 i=1,imt
          ztd(i,j) = ztd(i,j)*cpf(i,j)
790     continue
800   continue
c
c-----------------------------------------------------------------
c     estimate a first guess by extrapolating the two
c     previous solutions forward in time.
c-----------------------------------------------------------------
c
      fxa=c1
      if (mixts .or. mxpas2) fxa=p5
      if(.not.mxpas2)then
        do 900 j=1,jmt
        do 890 i=1,imt
          ptd(i,j,np) = fxa*(c2*ptd(i,j,nc)-ptd(i,j,nm))
890     continue
900     continue
c
c-----------------------------------------------------------------
c  allow for change in nm and nc indices during second pass
c  of euler backward
c-----------------------------------------------------------------
c
      else
        do 920 j=1,jmt
        do 910 i=1,imt
          ptd(i,j,np) = fxa*(c2*ptd(i,j,nm)-ptd(i,j,np))
910     continue
920     continue
      endif
c
c     compute criterion for convergence and set residuals to zero
c
      crtp = crit*fxa
      do 1000 j=1,jmt
      do 990 i=1,imt
        res(i,j) = c0
990   continue
1000  continue
c
c==================================================================
c     begin section to do the relaxation
c==================================================================
c
      nwcnt = 0
      nwmax = 300
      mscan = 0
1001  continue
      mscan = mscan + 1
      resmax = c0
c
c------------------------------------------------------------------
c     consider the arrays as being defined on the squares of an
c     "imt by jmt" checkerboard. take two passes: first solve the
c     equation on the black squares, then on the red squares.
c------------------------------------------------------------------
c
      do 1200 icolor=1,2
        do 1190 j=3,jscan
          ishypr   = 1 + mod(icolor+j,2)
          if (ishypr .eq. 1) ishypr = 3
          do 1160 i=ishypr,imtm1,2
            res(i,j) = rmsk(i,j)*(cfn(i,j)*ptd(i,j+1,np)
     &                 +cfs(i,j)*ptd(i,j-1,np)
     &                 +cfe(i,j)*ptd(i+1,j,np)
     &                 +cfw(i,j)*ptd(i-1,j,np)
     &                 -ptd(i,j,np)-ztd(i,j))
1160      continue
          if (ishypr .eq. 3) ishypr = 1
          res(1,j)   = res(imum1,j)
          res(imt,j) = res(2,j)
c
c------------------------------------------------------------------
c         make a correction to ptd based on the residuals
c------------------------------------------------------------------
c
          do 1170 i=ishypr,imt,2
            ptd(i,j,np) = ptd(i,j,np) + sor*res(i,j)
1170      continue
c
c------------------------------------------------------------------
c         find the maximum absolute residual to determine 
c         convergence
c------------------------------------------------------------------
c
          do 1180 i=ishypr,imt,2
            resmax = max(abs(res(i,j)),resmax)
1180      continue
1190    continue
1200  continue
c
c------------------------------------------------------------------
c     do hole relaxation for each island
c------------------------------------------------------------------
c
      do 1400 isle=1,nisle
        resis = c0
        do 1360 n=1,nippts(isle)
          i = iperm(iofs(isle)+n)
          j = jperm(iofs(isle)+n)
          resis = resis + (cfn(i,j)*ptd(i,j+1,np)
     &                  +  cfs(i,j)*ptd(i,j-1,np)
     &                  +  cfe(i,j)*ptd(i+1,j,np)
     &                  +  cfw(i,j)*ptd(i-1,j,np)
     &                  -  ptd(i,j,np) - ztd(i,j))*cof(i,j)
1360    continue
c
c       normalize the island residual and update the maximum
c       absolute residual of the relaxation if necessary
c
        resis  = resis*cofis(isle)
        resmax = max(abs(resis),resmax)
        resis  = resis*sor
c
c       make a correction to ptd over the island and its perimeter
c
        do 1370 n=1,nippts(isle)
          i = iperm(iofs(isle)+n)
          j = jperm(iofs(isle)+n)
          ptd(i,j,np) = ptd(i,j,np) + resis
1370    continue
1400  continue
c
c------------------------------------------------------------------
c  set cyclic boundary condition
c------------------------------------------------------------------
c
      do 1600 j=1,jmt
        ptd(1,j,np)   = ptd(imum1,j,np)
        ptd(imu,j,np) = ptd(2,j,np)
1600  continue
c
c------------------------------------------------------------------
c     test maximum residual for convergence of the relaxation.
c     if not converged, proceed with another scan.
c     (if the number of scans reaches mxscan, leave the loop)
c------------------------------------------------------------------
c
      if (resmax.ge.crtp .and. mscan.lt.mxscan) go to 1001
c
c==================================================================
c  end of the relaxation
c==================================================================
c
c------------------------------------------------------------------
c  update the stream function based on the relaxation solution
c------------------------------------------------------------------
c
      do 1700 j=1,jmt
      do 1690 i=1,imt
        p(i,j,np) = p(i,j,nm) + ptd(i,j,np)
1690  continue
1700  continue
c
c------------------------------------------------------------------
c    on a mixing timestep, alter ptd to be consistent with
c    normal, leap-frog stepping)
c------------------------------------------------------------------
c
      if (mixts .and. eb) return
      if (mxpas2 .or. mixts) then
        do 1900 j=1,jmt
        do 1890 i=1,imt
          ptd(i,j,np) = c2*ptd(i,j,np)
1890    continue
1900    continue
      endif
      return
      end

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
c
c------------------------------------------------------------------
c     set no flux condition for tracer at  bottom.
c     relax to "bcest" values with 50 day timescale at top
c------------------------------------------------------------------
c
      gamma = secday/50.0
      do 100 m=1,nt
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  DEBUG !! set flux to zero
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        stf(m)=c0
c        stf(m) = gamma*(stf(m)-t(1,ic,jc,m,nm))*dz(1)
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
c     references:
c       Hellerman, S, and M. Rosenstein, normal monthly wind stress
c         over the world ocean with error estimates, Journal of
c         Physical Oceanography, 13, 1093-1104,1983.
c       Levitus, S., Climatological atlas of the world ocean, NOAA
c         Prof. Paper 13, US Gov't printing Office, Washington, DC, 
c         1982.
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

      subroutine addv( i, j)
c
c==================================================================
c
c  add barotropic to baroclinic velocities
c  and calculate barotropic correction for
c      momentum advection equation
c
c==================================================================
c
#include "param.h"
#include "grdvar.h"
#include "emode.h"
#include "levind.h"
#include "timelv.h"
#include "slabs.h"
c
c-----------------------------------------------------------------
c  calculate barotropic velocity
c-----------------------------------------------------------------
c
      diag1  = p(i+1,j+1,np) - p(i  ,j  ,np)
      diag2  = p(i  ,j+1,np) - p(i+1,j  ,np)
      ubar   = -(diag1+diag2)*dy2r*hr(i,j)
      vbar   =  (diag1-diag2)*dx2r*hr(i,j)*csur(j)
c  add to baroclinic velocity
      do 100 k=1,kmu(i,j)
      q1temp=u(k,i,j,np)
      q2temp=v(k,i,j,np)
      u(k,i,j,np) = u(k,i,j,np)+ubar
      v(k,i,j,np) = v(k,i,j,np)+vbar
#ifdef debugprint
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  DEBEG
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(i.eq.60.and.j.eq.20.and.k.lt.4)then
      print *
      print *,' subroutine addv after adding ubar. k=',k
      print 7,' u(k,60,20,np),v(k,60,20,np) old =',
     &          q1temp,q2temp
      print 7,' ubar,vbar                       =',ubar,vbar
      print 7,' u(k,60,20,np),v(k,60,20,np) new =',
     &          u(k,60,20,np),v(k,60,20,np)
   7  format(a,1p8e16.6)
      endif
c!!END DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
 100  continue
c
      return
      end

      subroutine corv( i, j)
c
c==================================================================
c
c   calculate barotropic correction for momentum advection equation
c
c==================================================================
c
#include "param.h"
#include "grdvar.h"
#include "emode.h"
#include "levind.h"
#include "timelv.h"
#include "slabs.h"
c
c-----------------------------------------------------------------
c  calculate barotropic velocity
c-----------------------------------------------------------------
c
      diag1  = p(i+1,j+1,nc) - p(i  ,j  ,nc)
      diag2  = p(i  ,j+1,nc) - p(i+1,j  ,nc)
      ubar   = -(diag1+diag2)*dy2r*hr(i,j)
      vbar   =  (diag1-diag2)*dx2r*hr(i,j)*csur(j)
c  add to baroclinic velocity
c
c-----------------------------------------------------------------
c  calculate barotropic correction
c-----------------------------------------------------------------
c
c  calculate v barotropic velocity in box to south
      diag1  = p(i+1,j  ,nc) - p(i  ,j-1,nc)
      diag2  = p(i  ,j  ,nc) - p(i+1,j-1,nc)
      vbars  =  (diag1-diag2)*dx2r*hr(i,j-1)*csur(j-1)
c  calculate u barotropic velocity in box to west
      diag1  = p(i  ,j+1,nc) - p(i-1,j  ,nc)
      diag2  = p(i-1,j+1,nc) - p(i  ,j  ,nc)
      ubarw  = -(diag1+diag2)*dy2r*hr(i-1,j)
c  calculate external mode velocities
      sfu(i,j) = -(p(i,j+1,nc)-p(i,j,nc))*dyr*
     &                   min(hr(i-1,j),hr(i,j))
      sfv(i,j) =  (p(i+1,j,nc)-p(i,j,nc))*dxr*
     &                   min(hr(i,j),hr(i,j-1))*cstr(j)
c  calculate  correction
      q1temp = sfu(i,j)
      q2temp = sfv(i,j)
      sfu(i,j) = sfu(i,j) - p5*(ubar+ubarw)
      sfv(i,j) = sfv(i,j) - p5*(vbar+vbars)
c
#ifdef debugprint
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  DEBEG
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(i.eq.60.and.j.eq.20)then
      print *
      print *,' subroutine "corv" after calculating sfu,sfv'
      print 7,' ubar ,vbar    (60,20)     =',ubar,vbar
      print 7,' ubarw,vbars               =',ubarw,vbars
      print 7,' sfu(60,20),sfv(60,20) old =',q1temp,q2temp
      print 7,' sfu(60,20),sfv(60,20) new =',sfu(60,20),sfv(60,20)
   7  format(a,1p8e16.6)
      endif
c!!END DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
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
c  DEBUG !!
      data in/1/
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
C  DEBUG !!!! test m.eq.1 or 2
      if(m.eq.0.and.(t(k,ic,jc,m,np).ne.t(k,ic,jc,m,nm)))
     &  then
        in = in+1
        print *
        print *,'  diag k,ic,jc,m,np,nm=',k,ic,jc,m,np,nm
        print *,'  t(k,ic,jc,m,np),t(k,ic,jc,m,nm), ',
     &            't(k,ic,jc,m,np)-t(k,ic,jc,m,nm)'
      print *,     t(k,ic,jc,m,np),t(k,ic,jc,m,nm),
     &             t(k,ic,jc,m,np)-t(k,ic,jc,m,nm)
      print *,'  tddt(2),dtabs(2),tvar(2)='
      print *,tddt(2),dtabs(2),tvar(2)
        if(in.gt.100)stop
      endif
 500    continue
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
#include "emode.h"
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
c                                      to -2.0 and 40.0 over land.
c------------------------------------------------------------------
c
      do 200 n=1,3
      do 200 j=1,jmt
      do 200 i=1,imt
      do 200 k=1,km
        u(k,i,j,n)=c0
        v(k,i,j,n)=c0
        if(k.gt.kmt(i,j))then
          t(k,i,j,1,n)=-2.0
          t(k,i,j,2,n)=40.0
        else
          t(k,i,j,1,n)= theta0 (yt(j), zt(k), cst(j))
          t(k,i,j,2,n)= 34.9
        endif
200   continue
c
c------------------------------------------------------------------
c     initialize stream function p and rate of change ztd
c------------------------------------------------------------------
c
      do 300 n=1,3
        do 290 j=1,jmt
        do 280 i=1,imt
          p(i,j,n) = c0
          ptd(i,j,n) = c0
280     continue
290     continue
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
      subroutine iperim (kmz)
 
c
c==================================================================
c
c     find island perimeter coordinates & analyze for errors in 
c     land mass distribution
c
c     inputs:
c
c     kmz   = scratch area
c     kmt   = number of levels on the "t" grid
c     kmu   = number of levels on the "u" grid
c
c     output:
c
c     all island variables in "index.h"
c
c==================================================================
c
      character*1 dir, blank
      logical prblm, fatal
#include "param.h"
#include "coord.h"
#include "index.h"
#include "levind.h"
      dimension kmz(imt,jmt)
      dimension ncol(imt)
c
      data blank /' '/
c
c------------------------------------------------------------------
c     statement functions
c------------------------------------------------------------------
c
      iwest(i) = mod(i-2+imt,imt) + 1
      ieast(i) = mod(i,imt) + 1
c
c------------------------------------------------------------------
c     save the original "kmt" field (to be restored later)
c     and initialize "chisl" field
c------------------------------------------------------------------
c
      do 100 j=1,jmt
        do 90 i=1,imt
          kmz(i,j) = kmt(i,j)
90      continue
100   continue
c
c------------------------------------------------------------------
c     search for & knockout isolated bays ... "t" grid boxes at the
c     surface which can't be influenced by advection
c------------------------------------------------------------------
c
      do 400 j=2,jmtm1
        do 390 i=2,imtm1
          if (kmt(i,j) .ne. 0) then
            if (kmu(i,j) .eq. 0 .and. kmu(i-1,j) .eq. 0 .and.
     &        kmu(i,j-1) .eq. 0 .and. kmu(i-1,j-1) .eq. 0) then
              kmt(i,j) = 0
            endif
          endif
390     continue
400   continue
c
c------------------------------------------------------------------
c     find the coordinates of the perimeter of each island
c------------------------------------------------------------------
c
      fatal = .false.
      do 600 n=1,nisle
c
c------------------------------------------------------------------
c       convert island latitude & longitude "seed" values to nearest
c       model points incices.
c------------------------------------------------------------------
c
        i = indp (alonis(n), xt, imt)
        j = indp (alatis(n), yt, jmt)
        if (kmt(i,j) .ne. 0) then
          write (stdout,'(a17,i2,a18,5x,a5,f8.2,a6,f8.2)')
     &    ' => island seed #',n,' is an ocean point'
     &,   ' lon=',alonis(n),', lat=',alatis(n)
          stop '=>iperim'
        endif
c
c------------------------------------------------------------------
c       determine a first perimeter point
c------------------------------------------------------------------
c
        if (j .le. jmt/2) then
          jend = jmt
          jd   = 1
        else
          jend = 1
          jd   = -1
        endif
        do 500 jj=j,jend,jd
          if (kmt(i,jj) .ne. 0) then
            j = jj
            go to 510
          endif
500     continue
510     continue
        i0  = i
        j0  = j
        dir = blank
c
c------------------------------------------------------------------
c       scan for the rest of the perimeter
c------------------------------------------------------------------
c
        if (n .eq. 1) then
          iofs(n)  = 0
        else
          iofs(n)  = nippts(n-1) + iofs(n-1)
        endif
        lenp       = iofs(n) + 1
        iperm(lenp) = i
        jperm(lenp) = j
c
        do 580 m=1,maxipp
          iw  = iwest(i)
          ie  = ieast(i)
          js  = max(1,j-1)
          jn  = min(jmt,j+1)
          if (iw .eq. 1) then
            iw  = imt-1
          endif
          if (ie .eq. imt) then
            ie  = 2
          endif
c
c------------------------------------------------------------------
c         analyze the geometry to verify that the island and any other
c         land mass are separated by at least 2 ocean "t" grid boxes
c------------------------------------------------------------------
c
          prblm = .false.
c
          if ((kmt(ie,jn) .eq. 0 .and. kmt(iw,js) .eq. 0) .or.
     &        (kmt(ie,js) .eq. 0 .and. kmt(iw,jn) .eq. 0)) then
            if (kmt(iw,j) .ne. 0 .and. kmt(ie,j) .ne. 0 .and.
     &          kmt(i,jn) .ne. 0 .and. kmt(i,js) .ne. 0 ) prblm = .true.
            if (kmt(iw,j) .ne. 0 .and. kmu(iw,j) .eq. 0 .and.
     &          kmu(iw,js) .eq. 0) prblm = .true.
            if (kmt(ie,j) .ne. 0 .and. kmu(i,j) .eq. 0 .and.
     &          kmu(i,js) .eq. 0) prblm = .true.
            if (kmt(i,jn) .ne. 0 .and. kmu(i,j) .eq. 0 .and.
     &          kmu(iw,j) .eq. 0) prblm = .true.
            if (kmt(i,js) .ne. 0 .and. kmu(i,js) .eq. 0 .and.
     &          kmu(iw,js) .eq. 0) prblm = .true.
          endif
          if (prblm) then
            write (stdout,9200) n, i, j, xt(i), yt(j)
            fatal = .true.
          endif
c
c------------------------------------------------------------------
c         look west (only if not comming from west)
c------------------------------------------------------------------
c
          k2 = kmt(iw,jn)  * kmt(iw,js)
          k3 = kmt(i,jn)   * kmt(i,js)
          if (k2 .eq. 0 .or. k3 .eq. 0) then
            k = 0
          else
            k = 1
          endif
          if (kmt(iw,j) .ne. 0 .and. k .eq. 0 .and. dir .ne. 'e') then
            if (lenp .lt. maxipp) then
              if (iw .ne. i0 .or. j .ne. j0) then
                dir          = 'w'
                lenp         = lenp + 1
                i            = iw
                iperm(lenp)  = iw
                jperm(lenp)  = j
                go to 580
              else
                go to 590
              endif
            else
              write (stdout,'(a23)') '=> maxipp limit reached'
              stop '=>west'
            endif
          endif
c
c------------------------------------------------------------------
c         look north (only if not comming from north)
c------------------------------------------------------------------
c
          k2 = kmt(ie,jn)  * kmt(iw,jn)
          k3 = kmt(ie,j)   * kmt(iw,j)
          if (k2 .eq. 0 .or. k3 .eq. 0) then
            k = 0
          else
            k = 1
          endif
          if (kmt(i,jn) .ne. 0 .and. k .eq. 0 .and. jn .ne. j
     &        .and. dir .ne. 's') then
            if (lenp .lt. maxipp) then
              if (i .ne. i0 .or. jn .ne. j0) then
                dir          = 'n'
                lenp         = lenp + 1
                j            = jn
                iperm(lenp)  = i
                jperm(lenp)  = jn
                go to 580
              else
                go to 590
              endif
            else
              write (stdout,'(a23)') '=> maxipp limit reached'
              stop '=>north'
            endif
          endif
c
c------------------------------------------------------------------
c         look east (only if not comming from east)
c------------------------------------------------------------------
c
          k2 = kmt(ie,jn)  * kmt(ie,js)
          k3 = kmt(i,jn)   * kmt(i,js)
          if (k2 .eq. 0 .or. k3 .eq. 0) then
            k = 0
          else
            k = 1
          endif
          if (kmt(ie,j) .ne. 0 .and. k .eq. 0 .and. dir .ne. 'w') then
            if (lenp .lt. maxipp) then
              if (ie .ne. i0 .or. j .ne. j0) then
                dir          = 'e'
                lenp         = lenp + 1
                i            = ie
                iperm(lenp)  = ie
                jperm(lenp)  = j
                go to 580
              else
                go to 590
              endif
            else
              write (stdout,'(a23)') '=> maxipp limit reached'
              stop '=>east'
            endif
          endif
c
c------------------------------------------------------------------
c         look south (only if not comming from south)
c------------------------------------------------------------------
c
          k2 = kmt(ie,js)  * kmt(iw,js)
          k3 = kmt(ie,j)   * kmt(iw,j)
          if (k2 .eq. 0 .or. k3 .eq. 0) then
            k = 0
          else
            k = 1
          endif
          if (kmt(i,js) .ne. 0 .and. k .eq. 0 .and. js .ne. j
     &        .and. dir .ne. 'n') then
            if (lenp .lt. maxipp) then
              if (i .ne. i0 .or. js .ne. j0) then
                dir          = 's'
                lenp         = lenp + 1
                j            = js
                iperm(lenp)  = i
                jperm(lenp)  = js
                go to 580
              else
                go to 590
              endif
            else
              write (stdout,'(a23)') '=> maxipp limit reached'
              stop '=>south'
            endif
          endif
580     continue
590     continue
        nippts(n) = lenp - iofs(n)
        write (stdout,9000) n, nippts(n)
        ls = iofs(n) + 1
        le = lenp
        write (stdout,9100) (iperm(l),jperm(l), l=ls,le)
600   continue
c
c------------------------------------------------------------------
c     restore the original "kmt" field
c------------------------------------------------------------------
c
      do 700 j=1,jmt
        do 690 i=1,imt
          kmt(i,j) = kmz(i,j)
690     continue
700   continue
c
c------------------------------------------------------------------
c     if a fatal error was encountered print out error message and stop
c     otherwise load perimeter field into "kmz" and print it out
c------------------------------------------------------------------
c
      if (fatal) then
c
        write (stdout,*)
     &   ' Fatal error(s) in "iperim" land mass analysis'
        write (stdout,*)
     &    ' The "kmt" field must be fixed up before continuing ...'
        stop '=>iperim'
c
      else
c
c------------------------------------------------------------------
c     initialize scratch "kmz" field
c------------------------------------------------------------------
c
        do 800 j=1,jmt
          do 790 i=1,imt
            kmz(i,j) = 0
790       continue
800     continue
c
c------------------------------------------------------------------
c       identify perimeter points with island number
c------------------------------------------------------------------
c
        do 900 ni=1,nisle
          ls = iofs(ni) + 1
          le = ls + nippts(ni) - 1
          do 890 m=ls,le
            kmz(iperm(m),jperm(m)) = ni
890       continue
900     continue
c
c------------------------------------------------------------------
c       print out island perimeter field
c------------------------------------------------------------------
c
        linemx = 100
        linel  = 105
        line   = linel - 5
        if (line .gt. linemx) line = linemx
        nwr = (imt/line) + 1
c
        do 950 i=1,imt
          ncol(i) = mod(i,10)
950     continue
c
       write (stdout,9300)
        do 1000 n=1,nwr
          ia = 1 + (line*(n-1))
          ib = ia + line - 1
          if (ib .gt. imt) ib = imt
          write(stdout,9400) (ncol(i),i=ia,ib)
          do 990 jj=1,jmt
            jjj = jmt - jj + 1
            write(stdout,9500)jjj,(kmz(i,jjj),i=ia,ib)
990       continue
1000    continue
c
      endif
c
9000  format (/' Island # ',i2,' perimeter length is ',i4, ' points.'
     &, ' (i,j) coordinates are:'/)
9100  format (10(1x,'(',i4,',',i4,') '))
9200  format (///1x, ' ==> Error on island # ',i2,' at (i,j) = (',i4,','
     &,  i4,') ... (lon,lat) = (',f7.2,',',f7.2,')'/
     &, 1x,'    There must be at least two ocean points separating '
     &,    'an island and another land mass.'/)
9300  format (/,' Island perimeter locations '/)
9400  format('  i=>',100(i1))
9500  format(1x,i3,1x,100(i1))
c
      return
c
      end

      function theta0 (ydeg, depth, coslat)
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
      save coeft, tmin, tmax, reflat
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
c     coslat = abs(cos(pi*ydeg/180.))
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
c     itt    = itt + 1
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
      dtday = dt*secday
      last   = (totday .ge. days) 
      diagts = tmnset (dgnstc, totday, dtday)
      mixts  = (mod(itt,nmix) .eq. 1) 
      prntsi = tmnset (tsi, totday, dtday)
      snapts = tmnset (snaps, totday, dtday) .and. .not. last
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
#include "emode.h"
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
c-----------------------------------------------------------------------
c       print the stream function
c-----------------------------------------------------------------------
c
      scl=1.e12
      is = indp (slonxy, xt, imt)
      ie = indp (elonxy, xt, imt)
      js = indp (slatxy, yt, jmt)
      je = indp (elatxy, yt, jmt)
      write (stdout,8000) ' stream function (sverdrups)'
     &, itt, xt(is), xt(ie), yt(js), yt(je), scl
      call matrix (p(1,1,np), imt, is, ie, -js, -je, scl)
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
#include "iounit.h"
#include "slabs.h"
#include "timelv.h"
#include "emode.h"
      dimension t1(imt,jmt,km,nt), t2(imt,jmt,km,nt), 
     &          t3(imt,jmt,km,nt)
      dimension u1(imt,jmt,km),  u2(imt,jmt,km), u3(imt,jmt,km)
      dimension v1(imt,jmt,km),  v2(imt,jmt,km), v3(imt,jmt,km)
      dimension p1(imt,jmt),   p2(imt,jmt),   p3(imt,jmt)
      dimension ptd1(imt,jmt), ptd2(imt,jmt), ptd3(imt,jmt)
      character*80 file
      equivalence (t1(1,1,1,1),t(1,1,1,1,1))
      equivalence (t2(1,1,1,1),t(1,1,1,1,2))
      equivalence (t3(1,1,1,1),t(1,1,1,1,3))
      equivalence (u1(1,1,1),u(1,1,1,1))
      equivalence (u2(1,1,1),u(1,1,1,2))
      equivalence (u3(1,1,1),u(1,1,1,3))
      equivalence (v1(1,1,1),v(1,1,1,1))
      equivalence (v2(1,1,1),v(1,1,1,2))
      equivalence (v3(1,1,1),v(1,1,1,3))
      equivalence (p1(1,1),p(1,1,1))
      equivalence (p2(1,1),p(1,1,2))
      equivalence (p3(1,1),p(1,1,3))
      equivalence (ptd1(1,1),ptd(1,1,1))
      equivalence (ptd2(1,1),ptd(1,1,2))
      equivalence (ptd3(1,1),ptd(1,1,3))
c
c     read restart data set from archive for timesteps itt 
c
      entry rdrest
c  check np equals 3
      if(np.ne.3)then
        write(stdout,21) np, nc, nm
        stop
      endif
c
c  read archived data.  Use file fnrest if defined.
c
      if(fnrest(1:1).ne.' ')then
      open (iorest, access='SEQUENTIAL', form='UNFORMATTED')
        open(iorest, file=fnrest, form='UNFORMATTED', 
     &       access='SEQUENTIAL')
      else
        open(iorest, form='UNFORMATTED', access='SEQUENTIAL')
      endif
      read(iorest) itt0, totsec, totday, years, kmt,
     &             p3, ptd3, ptd2, t3, u3, v3
      close(iorest)
      write (stdout,22) iorest, itt0
      if(itt.ne.itt0)then
        write(stdout,23)itt
        stop
      endif
      return
c
  21  format(' Routine restio.',/,
     &       ' Trying to restart with (np not equal to 3).',/,
     &       ' np, nc, nm  =',3i5,/,' Programme stopped.')
  22  format(' Routine restio',/,
     &       ' Restart file read from unit ',i2,'.  timestep =',i7)
  23  format(' Routine restio',/,
     &       ' Timestep of restart file not equal to requested',
     &       ' timestep.  Requested timestep =',i7,/,
     &       ' Program stopped')
c
c     write restart data set to archive.
c
      entry wrrest
c
c  create archive file name
c
      write (file,'(i10.10)')nint(totday)
      file = 'd'//file(7:10)//'.restart'
      print *,'**file**',file,'*****'
      open(iorest, file=file, form='UNFORMATTED',
     &       access='SEQUENTIAL')
c
      write (stdout,24) itt, file, iorest
      write(iorest) itt, totsec, totday, years, kmt
      if(np.eq.1) write(iorest)  p1, ptd1
      if(np.eq.2) write(iorest)  p2, ptd2
      if(np.eq.3) write(iorest)  p3, ptd3
c
      if(nc.eq.1) write(iorest)  ptd1
      if(nc.eq.2) write(iorest)  ptd2
      if(nc.eq.3) write(iorest)  ptd3
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

