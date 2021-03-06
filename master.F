#ifndef Masterslave
      program master
#else
      subroutine master
#endif
c
c=======================================================================
c @(#) SCCS module: master.F  version: 1.1
c      Creation date: 10/13/97
c=======================================================================
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
c           David J. Webb          djw@soc.soton.ac.uk
c           Beverly A. de Cuevas   bac@soc.soton.ac.uk
c           Andrew C. Coward       acc@soc.soton.ac.uk
c           Catherine S. Gwilliam  cgw@soc.soton.ac.uk
c
c  Advice on CRAY T3D specifics has also been received from:
c
c           Mike O'Neill        meon@cray.com
c           Bob Carruthers      Bob.Carruthers@cray.com
c
c  Copyright 1993, 1994, 1995, 1996 and 1997.  
c                     D.J. Webb, B.A. de Cuevas, A.C. Coward, 
c                     A.C. Coward and C.S. Richmond,
c                     Southampton Oceanography Centre,
c                     Empress Dock, Southampton SO14 3ZH, U.K..
c
c  The code may be freely adapted and used without charge for
c  non-commercial purposes.  Publications which report on work that
c  has made use of the code should reference one or more of the
c  relevant papers produced by the OCCAM core team.
c  
c=======================================================================
c
#include "def_master.h"
#include "param.h"
#include "mesdta.h"
#include "scalar.h"
#include "switch.h"
#include "iounit.h"
#include "chmix.h"
#include "cvbc.h"
#include "cvmix.h"
#include "coord.h"
#include "ctmngr.h"
#include "grdvar.h"
#include "timelv.h"
#include "levind.h"
#include "archive.h"
#include "varinfo.h"
#include "snaps.h"
#include "versno.h"
#ifndef Masterslave
c
c  use MPI timing routines
c
      double precision total_time
#endif
c
      logical log1
      character*12 nodename
      character*8 archit
c
#include "moddata.h"
c
#ifdef no_namelist
      character*80 invar
#else
      namelist /contrl/ init, fnrest, days, restrt, nmix, eb
     &,         ncon, tsi,  dgnstc, snapd, archd, acor, idebug
# ifdef de_checkbd
     &,                 dchkbd
# endif
     &,         ftrest, ftsnap, ftarch, ispvar, isplev
     &,         statms, statis
      namelist /tsteps/ dtts, dtuv, dtbt
      namelist /eddy/   am, ah, fkpm, fkph, cdbot
#endif

c
c=======================================================================
c     begin
c     if necessary initialise MPI, start timing  
c      and initialise message passing logic
c     open input file and write version information
c=======================================================================
c
#ifndef Masterslave
      call MPI_INIT(ierr)
      start_time = MPI_WTIME()
      call MPI_BUFFER_ATTACH(bufmpi,LENBUF*NUMBUF*8,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,me,ierr)
#endif
c
      model(1) =' Stripped down MOM code for array processor.'
      model(2) =' Module: master.F.  Version: 1 '//
     &          ' Date: 10/13/97'
#ifdef no_namelist
      open (stdin,file='ocean.in2')
#else
      open (stdin,file='ocean.in')
#endif
      write (stdout,'(/2(25x,a80/)/)') model
c
c------------------------------------------------------------------
c     initialise master parameters
c------------------------------------------------------------------
c
      stlon = 0.0
      stlat = -72.0
      dxdeg = 4.0
      dydeg = 4.0
c
      init = .true.
      first= .true.
      days  = 2.0
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
      ftrest = 'hdf'
      ftsnap = 'hdf'
      ftarch = 'hdf'
      iorest = 21
      iosnap = 22
      statms = 10
      statis = 10
      idebug = 3
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
#ifdef de_checkbd
      dchkbd = 0.1
#endif
c
      pi     = c4*atan(c1)
      radian = c360/(c2*pi)
      omega  = pi/43082.0
      grav =980.6
      radius =6370.e5
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
c------------------------------------------------------------------
c     read in run parameters
c------------------------------------------------------------------
c
#ifdef no_namelist
      read (stdin,'(////)')
      read (stdin,*) invar	
      read (stdin,*) init, fnrest
c
      read (stdin,'( )')
      read (stdin,*) invar
      read (stdin,*) days, restrt
c
      read (stdin,'( )') 
      read (stdin,*) invar
      read (stdin,*) nmix, eb, ncon, acor
c
      read (stdin,'( )')
      read (stdin,*) invar
      read (stdin,*) tsi, snapd, archd, idebug
c
      read (stdin,'( )')
      read (stdin,*) invar
      read (stdin,*) ftrest, ftsnap, ftarch
c
      read (stdin,'( )')
      read (stdin,*) invar
      read (stdin,*) dgnstc, dchkbd
c
      read (stdin,'( )')
      read (stdin,*) invar
      read (stdin,*) statms, statis
c
      read (stdin,'( )')
      read (stdin,*) invar
      read (stdin,*) am, ah, fkpm, fkph, cdbot
c
      read (stdin,'( )')
      read (stdin,*) invar
      read (stdin,*) dtts, dtuv, dtbt
c
      write(stdout,1450) init,fnrest(1:lnblnk(ftrest)),days,restrt,
     &                   nmix,eb,ncon,acor,tsi,snapd,archd,idebug,
     &                   ftrest,ftsnap,ftarch,
     &                   dgnstc,dchkbd,statms,statis
      write(stdout,1451) am,ah,fkpm,fkph,cdbot
      write(stdout,1452) dtts,dtuv,dtbt
# ifdef de_checkbd
      write(stdout,1453) dchkbd
1453  format(' dchkbd = ',g10.4 )
# endif
c
1450  format(' init  = ',l1,' fnrest= ',a,/,
     &       ' days  = ',g10.4,' restrt= ',l1,/,
     &       ' nmix  = ',i4,' eb    = ',l1,' ncon  = ',i4,
     &       ' acor  = ',g10.4,/,
     &       ' tsi   = ',g10.4,' snapd = ',g10.4,' archd = ',g10.4,
     &       ' idebug= ',i10,/,
     &       ' ftrest = ',a,' ftsnap = ',a,' ftarch = ',a,/,
     &       ' dgnstc = ',g10.4,' dchkbd = ',g10.4,/,
     &       ' statms = ',g10.4,' statis = ',g10.4)
1451  format(' am    = ',g10.4,' ah    = ',g10.4,/,' fkpm  = ',g10.4,
     &       ' fkph  = ',g10.4,/,' cdbot = ',g10.4)
1452  format(' dtts  = ',g10.4,' dtuv  = ',g10.4,' dtbt  = ',g10.4)
#else
      read  (stdin, contrl)
      write (stdout,contrl)
      read  (stdin, eddy)
      write (stdout, eddy)
      read  (stdin, tsteps)
      write (stdout,tsteps)
#endif
      close (unit=stdin)
c
c-----------------------------------------------------------------------
c     calculate free surface model (barotropic) timestep
c-----------------------------------------------------------------------
c
      ntbt = max(1,nint(dtuv/dtbt))
      if(dtuv.ne.dtbt*ntbt)then
        write(stdout,'(a)') ' dtuv not integer multiple of dtbt.'
        write(stdout,*) ' dtuv is  ',dtuv
        write(stdout,*) ' dtbt was ',dtbt
        dtbt = dtuv/ntbt
        write(stdout,*) ' dtbt now set to set to:',dtbt
      endif
c
c------------------------------------------------------------------
c     list code options chosen by user for this model run
c------------------------------------------------------------------
c
      call docmnt (stdout, .false., .false., .false., .false.)
c
c-----------------------------------------------------------------------
c     do all consistency checks now
c-----------------------------------------------------------------------
c
      call checks
c
c-----------------------------------------------------------------------
c    initialise grid arrays
c    initialise depths, kmt, and associated arrays
c    initialise processor map, kpn, associated arrays and nproc
c-----------------------------------------------------------------------
c
      call m_grids
      call init_kmt
      call init_kpn
#ifdef de_checkbd
c
c scale the checkerboard filter weighting factor
c
      dchkbd = dchkbd*grav*dtbt*zw(KM)
#endif
c
c-----------------------------------------------------------------------
c    Check MPI environment
c-----------------------------------------------------------------------
c
      call MPI_COMM_GROUP(MPI_COMM_WORLD,mpigroup,ierr)
      call MPI_GROUP_SIZE(mpigroup,iproc,ierr)
      if(iproc.ge.nproc+1)then
        write(stdout,101)nproc
      else
        write(stdout,102)nproc+1,iproc
        call MPI_ABORT(MPI_COMM_WORLD,-1)
      endif
c
      write(stdout,103)
 101  format(' Starting ocean model on ',i4,' processors')
 102  format(' Number processors required     = ',i5,/,
     &       ' Number of processors available = ',i5,/,
     &       ' Program stopping ... ',/,)
 103  format(' ')
c
c-----------------------------------------------------------------------
c  initialise message info arrays and broadcast salve ids.
c-----------------------------------------------------------------------
c
      call m_initmess
      call m_send(MSG_IDS)
c
c-----------------------------------------------------------------------
c broadcast the control/namelist and grid information
c-----------------------------------------------------------------------
c
      call m_send(MSG_CONTROL)
c
c-----------------------------------------------------------------------
c  broadcast array lengths and topographic information
c-----------------------------------------------------------------------
c
      call m_send(MSG_TOPOG)
c
c-----------------------------------------------------------------------
c initialise prininfo
c-----------------------------------------------------------------------
c
      call prininfo(1)
      call prininfo(3)
c
c-----------------------------------------------------------------------
c  read in and broadcast the restart data.  the slaves checkpoint, 
c  receiving the data, until they receive MSG_CONT
c-----------------------------------------------------------------------
c
      if(.not.init) then
        call m_send(MSG_TIMEVAR)
        call archrd(1)
        call m_send(MSG_CONT)
      endif
c
c-----------------------------------------------------------------------
c  initialise and broadcast the flux data sets
c-----------------------------------------------------------------------
c
      call initmet
      nyears = totday/365.25
      tday = totday - 365.25*nyears
      do 20 n=0,12
        if(tday.ge.eday(n).and.tday.lt.eday(n+1))then
          mlast = n
          goto 21
        endif
  20  continue
      print *,' error setting up met fields. tday = ',tday
      print *,' eday = ',eday
      call m_abort
c
  21  if(mlast.eq.0) mlast = 12
      call metrd(mlast)
      call m_send(MSG_MET1)
c
      mnext = mod(mlast,12) + 1
      call metrd(mnext)
      call m_send(MSG_MET2)
c
      mnextp = mod(mnext,12) + 1
      call metrd(mnextp)
      lmetp = .true.
      lmetq = .true.
      nmeta = nproc
      call m_send(MSG_MET3)
c
      call mchkpnt(7)
c
c-----------------------------------------------------------------------
c  broadcast the Levitus data sets
c-----------------------------------------------------------------------
c

c
c-----------------------------------------------------------------------
c  checkpoint - wait for all slaves to check in and then start.
c-----------------------------------------------------------------------
c
      call prininfo(1)
      call prininfo(3)
      call mchkpnt(10)
c
c-----------------------------------------------------------------------
c  enter timestepping loop
c-----------------------------------------------------------------------
c
      do 3400 loop=1,9999999
c
c-----------------------------------------------------------------------
c  checkpoint at beginning of timestep
c-----------------------------------------------------------------------
c
#ifdef SYNC_3
      call mchkpnt(11)
#endif
c
c-----------------------------------------------------------------------
c       update timestep, set time dependent logical switches
c       to determine program flow for timestep itt, itt-1 and itt-2.
c       note: timestep itt and time refer to the time of the new
c       fields being calculated!
c-----------------------------------------------------------------------
c
      first = loop.eq.1
      call tmngr (dtts)
c
c-----------------------------------------------------------------------
c  if flags set start processing
c-----------------------------------------------------------------------
c
      if(prntsi) call inittsi
      if(snapts) call initsnap
      if(archts) call initarch
      if(metts)  call nextmet
c
c-----------------------------------------------------------------------
c  if an archive is in progress and if previous, partial archive request 
c  has been fulfilled, then send off the next part:
c-----------------------------------------------------------------------
c
      if(larchp.and..not.larchq) call archreq
c
c-----------------------------------------------------------------------
c  if a snapshot is in progress and if previous, snapshot request 
c  has been fulfilled, then send off the next part:
c-----------------------------------------------------------------------
c
      if(lsnapp.and..not.lsnapq) call snapreq
c
c-----------------------------------------------------------------------
c  loop processing messages until safe to proceed
c  print information arrays every few seconds while waiting
c-----------------------------------------------------------------------
c
      call prininfo(3)
      call m_recv(MSG_ANY)
c +->
c |
  25  if((prntsn .and. ltsip) .or. (snaptn .and. lsnapp)
     &                        .or. (archtn .and. larchp)
     &                        .or. (lmettn .and. lmetp)
     &                        .or. (last   .and. ltsip )
     &                        .or. (last   .and. lsnapp)
     &                        .or. (last   .and. larchp) )then
#ifndef cray-t3d
        call sleep(1)
#endif
        call m_recv(MSG_ANY)
        call prininfo(2)
c |
c-|----------------------------------------------------------------
c | if an archive is in progress and if previous, partial archive 
c | request has been fulfilled, then send off the next part:
c-|----------------------------------------------------------------
c |
        if(larchp.and..not.larchq) call archreq
c |
c-|----------------------------------------------------------------
c | if a snapshot is in progress and if previous snapshot request 
c | has been fulfilled, then send off the next part:
c-|----------------------------------------------------------------
c |
        if(lsnapp.and..not.lsnapq)  call snapreq
c |
        goto 25
      endif      
c |
c +->
c
c-----------------------------------------------------------------------
c  end timestepping loop
c-----------------------------------------------------------------------
c
      if(last)goto 3401
      call prininfo(2)
3400  continue
3401  continue
c 
c     wait for slaves to reach end of run then stop
c
  30  call mchkpnt(40)
c
c-----------------------------------------------------------------------
c     End user program 
c-----------------------------------------------------------------------
c     documentation info at end of run
c-----------------------------------------------------------------------
c
      call docmnt (stdout, .true., .true., .true., .true.)
c
c-----------------------------------------------------------------------
c     leave MPI before exiting 
c-----------------------------------------------------------------------
c
#ifdef Masterslave
      return
#else 
c
c  write out total time and stop
c
      total_time = MPI_WTIME() - start_time
      write( stdout,9603)total_time
      write (stdout,9604)
      call MPI_FINALIZE(ierr)
      stop
9603  format('Total elapsed time for master is ',f10.4)
9604  format(1x,' ==> END of model run.',//)
#endif
      end

