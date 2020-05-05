#ifndef Masterslave
      program slave
#else
      subroutine slave
#endif
c
c @(#) SCCS module: slave.F  version: 1.1
c     Creation date: 03/16/95
c
c======================================================================
c
c  This program and the corresponding master program make up  
c  'moma.pvm', the pvm message passing versions of David Webb's 
c  moma code, itself based on the GFDL modular oceam model (MOM).
c  Master and slave have been developed as part of the UK's Ocean 
c  Circulation & Climate Advanced Modelling project (OCCAM).
c
c  The master program spawns the slave processes and then acts as
c  a supervisor and input-output manager.  It collects and collates
c  timestep information, archive data and status reports.  The 
c  master also responds to the slave's requests for met data, 
c  initial conditions and set-up information. 
c 
c  Each slave is responsible for its own communications with 
c  its neighbours and responds to the masters's requests for 
c  snapshot and archive data.
c 
c  The master and slave programs have been developed by members of 
c  the OCCAM core team including:
c
c           David J. Webb       djw@ios.ac.uk 
c           Beverly de Cuevas   bac@ios.ac.uk 
c           Andrew C. Coward    acc@ios.ac.uk 
c           Catherine Gwilliam  cgw@ios.ac.uk 
c 
c  Advice on CRAY T3D specifics has also been received from:
c
c           Mike O'Neill        meon@cray.com
c
c
c  Copyright 1993, 1994, 1995.  D.J. Webb.  
c                   Institute of Oceanographic Sciences,
c                   Godalming, Surrey GU8 5UB, U.K..
c
c  The code may be freely adapted and used without charge
c  for non-commercial purposes.
c  
c======================================================================
c
c  Slave consists of a message-passing wrapper around and within 
c  David Webb's MOMA program, Version 1.12.  Moma is a free 
c  surface ocean model designed for use with array processors.
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
c      pressure field is included using flag 'presetp'.
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
c=======================================================================
c
#include "def_slave.h"
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
#include "cvmix.h"
#include "frees.h"
#include "iounit.h"
#include "varinfo.h"
#include "versno.h"
#include "mesdta.h"
c
#ifndef Masterslave
# include "moddata.h"
#endif 
c
      character fnnum*3,fnslave*80
c
#ifndef no_namelist
      namelist /contrl/ init, fnrest, days, restrt, nmix, eb, ncon
     &,                 tsi,  dgnstc, snapd, archd, acor, idebug
      namelist /tsteps/ dtts, dtuv, dtbt
      namelist /eddy/   am, ah, fkpm, fkph, cdbot
# ifdef de_checkbd
     &,                 dchkbd
# endif
#endif
c
c
c=======================================================================
c
c  start timing and initialise message passing logic
c 
c=======================================================================
c
#ifndef Masterslave
      utime(3)=dtime(utime)
#endif
      call initmess
c
c     constants usually initialised in block data
c
#ifndef Masterslave
      call getmyid(mytid)
      call getmasterid(mtid)
#endif
c
c receive tid array from master and decide which slave we are
c
      call s_recv(MSG_IDS)
c
c     Determine which slave I am (1 -- nproc)
c
#ifdef cray-t3d
      if(nproc .gt. n$pes-1) then
        write(stdout,*)'nproc is larger than the number of PEs ',
     &                'available'
        write(stdout,*)'slave ', me,' stopping ...'
        call closedown(info)
        return
      endif
#endif
      do 57 i=1, nproc
         if( itida(i) .eq. mytid ) me = i
 57   continue
#ifdef cray-t3d
      if(me .gt. nproc) then
        write(stdout,*)'PE no ',me,' closing down - not needed '
        write(stdout,*)'PE no ',me,' nproc = ',nproc
        call closedown(info,me)
        return
      endif
#endif
c
c=======================================================================
c     begin introductory section
c     open output file and write MOM version information
c=======================================================================
c
c
c  create slave output file 
c
      ioslave = 23
      if(ioslave.ne.stdout)then
        write(fnnum,'(i3.3)') me
        fnslave = 'slave'//fnnum
        open (unit=ioslave,status='UNKNOWN',file=fnslave,
     &         iostat = iostat)
        if(iostat.ne.0)then
          write(stdout,*)' Error opening slave output file:',
     &                     fnslave
          write(stdout,*)' iostat =',iostat
          write(stdout,*)' slave using stdout ...'
          ioslave = stdout
        else
          rewind(ioslave)
        endif
      endif
c
      model(1) =' Stripped down MOM code for array processor.'
      model(2) =' Module: slave.F.  Version: ????. '//
     &          ' Date:       16:10'
      print *, model
      print *,ioslave, stdout
      write(ioslave,*)model
      write(ioslave,*)'  me, mtid, mytid = ',me,mtid,mytid
      write(ioslave,*)'  itida =',(itida(n),n=1,nproc)
c
c=======================================================================
c   constants
c=======================================================================
c
      pi     = c4*atan(c1)
      radian = c360/(c2*pi)
      omega  = pi/43082.0
      grav =980.6
      radius =6370.e5
c
c=======================================================================
c   receive basic control information
c=======================================================================
c
      call s_recv(MSG_CONTROL)
#ifdef no_namelist
      write(ioslave,1452) init,fnrest,days,restrt,nmix,eb,ncon,
     &                   tsi,dgnstc,snapd,archd,acor,idebug
      write(ioslave,1451) dtts,dtuv,dtbt
      write(ioslave,1450) am,ah,fkpm,fkph,cdbot
#ifdef de_checkbd
      write(stdout,1448) dchkbd
1448  format(' dchkbd = ',g10.4 )
#endif
c
1450  format(' am    = ',g10.4,' ah    = ',g10.4 / ' fkpm  = ',g10.4,
     &       ' fkph  = ',g10.4 / ' cdbot = ',g10.4)
1451  format(' dtts  = ',g10.4,' dtuv  = ',g10.4,' dtbt  = ',g10.4)
1452  format(' init  = ',l1,' fnrest= ',a / ' days  = ',g10.4,
     &       ' restrt= ',l1,' nmix  = ',i4,' eb    = ',l1,' ncon  = ',i4
     &     / ' tsi   = ',g10.4,' dgnstc= ',g10.4,' snapd = ',g10.4
     &     / ' archd = ',g10.4,' acor  = ',g10.4,' idebug= ',i10)
#else
      write (ioslave,contrl)
      write (ioslave, eddy)
      write (ioslave,tsteps)
#endif
c
c
c=======================================================================
c     end introductory section
c=======================================================================
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
c=======================================================================
c     receive the number of vertical levels on the "t" grid
c     and grid related information
c=======================================================================
c
      call s_recv(MSG_TOPOG)
      itm1 = it-1
      jtm1 = jt-1
c
c-----------------------------------------------------------------------
c       print map of "kmt" levels.
c-----------------------------------------------------------------------
c
      if(.true.)then
        write(ioslave,*) 'Kmt map>>>>>>>>>>>>>>>>>>>'
        do 1300 ibk=1,it,40
          isp  = ibk
          iept = ibk + 40 - 1
          if(iept.gt.it) iept=it
          write (ioslave,'(/, 4x, 40i3)') (ii, ii=isp,iept)
          do 1290 jrev=1,jt
            j=jt-jrev+1
            write (ioslave,'(1x,i3, 40i3)')j,(kmt(i,j),i=isp,iept)
1290      continue
1300    continue
      endif
c
c-----------------------------------------------------------------------
c       print map of "kmu" levels. 
c-----------------------------------------------------------------------
c 
      if(.false.)then
        write(ioslave,*) 'Kmu map>>>>>>>>>>>>>>>>>>>'
        do 1301 ibk=1,it,40 
          isp  = ibk   
          iept = ibk + 40 - 1 
          if(iept.gt.it) iept=it
          write (ioslave,'(/, 4x, 40i3)') (ii, ii=isp,iept) 
          do 1291 jrev=1,jt
            j=jt-jrev+1  
            write (ioslave,'(1x,i3, 40i3)')j,(kmu(i,j),i=isp,iept) 
1291      continue 
1301    continue 
      endif
c
c-----------------------------------------------------------------------
c       print map of "kpn" array
c-----------------------------------------------------------------------
c
      if(.false.)then
        write(ioslave,*) 'kpn map>>>>>>>>>>>>>>>>>>>'
        do 1320 ibk=1,it,40
          isp  = ibk
          iept = ibk + 40 - 1
          if(iept.gt.it) iept=it
          write (ioslave,'(/, 4x, 40i3)') (ii, ii=isp,iept)
          do 1315 jrev=1,jt
            j=jt-jrev+1
            write (ioslave,'(1x,i3, 40i3)')j,(kpn(i,j),i=isp,iept)
1315      continue
1320    continue
      endif
c
      call bound
      call s_grids
c
c-----------------------------------------------------------------------
c     is this a start from initial conditions?
c-----------------------------------------------------------------------
c
      if (init) then
        call ocn1st
      else
c
c-----------------------------------------------------------------------
c checkpoint and process messages until the master broadcasts a
c continue. The master will keep the slaves at the checkpoint until
c all the restart data has been received and stored in the buffer
c-----------------------------------------------------------------------
c
        call schkpnt(4)
c
c transfer data from the buffer storage
c
        call rdrest
      endif
c
c-----------------------------------------------------------------------
c     compute maximum number of vertical levels on the "u" grid
c-----------------------------------------------------------------------
c
#ifdef oldadv 
        do 1110 j=2,jtm1
        do 1110 i=2,itm1
          kmd(i,j)=max(kmu(i+1,j),kmu(i-1,j),kmu(i,j+1), 
     &                  kmu(i,j-1),kmu(i,j)) 
1110    continue   
#else 
        do 1120 j=2,jtm1
        do 1120 i=2,itm1
          kmd(i,j)=max(kmu(i-1,j-1),kmu(i,j-1),kmu(i+1,j-1),
     &                  kmu(i-1,j),kmu(i,j),kmu(i+1,j),
     &                  kmu(i-1,j+1),kmu(i,j+1),kmu(i+1,j+1))
1120    continue
#endif
c
c-----------------------------------------------------------------------
c     compute area and volume of ocean ("t,s" grid boxes)
c-----------------------------------------------------------------------
c
      area   = c0
      volume = c0
      ocnp   = 0
c
      do 700 j=2,jtm1
      do 690 i=2,itm1
      if (kmt(i,j) .gt. 0 .and. kpn(i,j).eq.me) then
        area   = area + cst(j)*dx*dy
        volume = volume + cst(j)*dx*dy*zw(kmt(i,j))
        ocnp   = ocnp + float(kmt(i,j))
      endif
690   continue
700   continue
      write (ioslave,9341) area, volume
c
c-----------------------------------------------------------------------
c     compute depths and reciprocal depths
c-----------------------------------------------------------------------
c
      do 1400 j=1,jt
      do 1390 i=1,it
        h(i,j)  = c0
        hr(i,j) = c0
        if (kmu(i,j) .ne. 0) then
          h (i,j) = zw(kmu(i,j))
          hr(i,j) = c1/zw(kmu(i,j))
        endif
1390  continue
1400  continue
c
c-----------------------------------------------------------------------
c     initialize various things
c-----------------------------------------------------------------------
c
      do 2600 j=1,jt
      do 2590 i=1,it
        zu(i,j)  = c0
        zv(i,j)  = c0
2590  continue
2600  continue
c
c     initialize the coriolis factor
c
      do 2660 j=1,jt
        fcor(j) = c2*omega*sine(j)
2660  continue
c
c     initialize the diffusion factors
c
      do 2670 j=2,jtm1
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
c  report to master and wait
c
      call schkpnt(10)
c
c-----------------------------------------------------------------------
c     start the time step loop
c-----------------------------------------------------------------------
c
      do 3400 loop=1,9999999
c
c  call checkpoint and clear input buffers
c
#ifdef SYNC_3
        call schkpnt(11)
#else
        call sendposn(11)
#endif
        call s_recv(MSG_ANY)
c
c-----------------------------------------------------------------------
c       update timestep, set time dependent logical switches
c       to determine program flow for timestep itt, itt-1 and itt-2.
c       note: timestep itt and time refer to the time of the new
c       fields being calculated!
c-----------------------------------------------------------------------
c
        first= loop .eq. 1
        call tmngr (dtts)
        call step 
c
c-----------------------------------------------------------------------
c      call routine to check requests from master, then send 
c      time step information, snapshots and archive data 
c      as required.  Loop until safe to continue
c-----------------------------------------------------------------------
c
        call sendposn(35)
 25   continue
      call s_recv(MSG_ANY)
      call flush(ioslave)
      if(prntsn .and. ltsip)then
        call sleep(1)
        goto 25
      endif
      if(archtn .and. larchp)then
c
c Pause to allow master request to come in
c
#ifndef cray-t3d
       if(.not.larchq) call sleep(1)
#endif
       goto 25
      endif
      if(snaptn .and. lsnapp)then
c
c Pause to allow master request to come in
c
#ifndef cray-t3d
       if(.not.lsnapq) call sleep(1)
#endif
       goto 25 
      endif
c
#ifdef cray-t3d
      time_now=rtc()
      utime(1)=time_now-utime(2)
      utime(2) = time_now
      time_elapsed =(time_now-utime(3))/1.5e8
      if (prntsi) write (ioslave,9602) itt, stamp, ektot/volume,
     &        dtabs(1)/volume, dtabs(2)/volume,time_elapsed
#else
      utime(3)=dtime(utime)+utime(3)
      if (prntsi) write (ioslave,9602) itt, stamp, ektot/volume,
     &        dtabs(1)/volume, dtabs(2)/volume,utime(1),utime(3)
#endif
      call flush(ioslave)
c
      if (diagts) call pdiag
      if (archts) then
       call wrrest
       larchp = .true.
      endif
      if (snapts) then
       call snap
       lsnapp = .true.
      endif
c
c------------------------------------------------------------------
c       end timestepping loop
c------------------------------------------------------------------
c
      if (last) go to 3401
3400  continue
3401  continue
c
c-----------------------------------------------------------------
c checkpoint until all outstanding requests have been processed.
c The checkpoint routine will spin receiving and acting upon
c messages until told to proceed by caesar. Additionally, s_recv
c (called by schkpnt) will send out timestep, archive and snapshot
c data still awaiting collection by caesar.
c-----------------------------------------------------------------
      call schkpnt(40)
c
c------------------------------------------------------------------
c     close all units  (stderr is not closed because it is set equal
c     to ioslave in pconst.h).
c------------------------------------------------------------------
c
      write (ioslave,9603) 
      close (unit=ioslave,status='keep')
      if(ioslave.ne.stdout)then
        close (unit=stdout,status='keep')
      endif
      call closedown(info)
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
c     &        ' dtemp=',1pe13.6,' dsalt=',1pe13.6)
     &        ' dtemp=',1pe13.6,' dsalt=',1pe13.6,' times=',0p2f10.2)
9603  format(1x,' ==> END of model run.',//)
      end

