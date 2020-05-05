      subroutine m_initmess
c
c=======================================================================
c @(#) SCCS module: m_msg.F  version: 1.8
c     Creation date: 08/23/95
c=======================================================================
c
c
c  subroutine used in master to:
c
c  1.  initialsie info arrays for message passing.
c  2.  initialise flags and logical variables for
c      printing timestep information, and saving 
c      snapshot and archive information.
c
c=======================================================================
c
#include "def_master.h"
#include "param.h"
#include "levind.h"
#include "mesdta.h"
c
c=======================================================================
c  1.  initialise message passing arrays
c=======================================================================
c
      do 10 j=1,MXSLAVE
      do 10 i=1,3
        infoa(i,j)=0
  10  continue
c
c=======================================================================
c  2.  initialise tsi, snapshot, archive, and flux flags
c=======================================================================
c
      do 20 i=1,MXSLAVE
        ltsia(i)  = .false.
        lsnapa(i) = .false.
        lfluxa(i) = .true.
  20  continue
      ltsip  = .false.
      lsnapp = .false.
      lsnapq = .false.
      larchp = .false.
      larchq = .false.
c
      return
      end

      subroutine mchkpnt(iposn0)
c
c=======================================================================
c  master subroutine to checkpoint until if all slaves have  
c  reached position iposn on current timestep and then 
c  issue message to continue.
c
c  input:  iposn - position for checking
c=======================================================================
c
#include "param.h"
#include "ctmngr.h"
#include "mesdta.h"
c
c  set position
c
      iposn = iposn0
c
c  receive and process all waiting messages
c
  10  call m_recv(MSG_ANY) 
c
c  check position and timestep correct for all processors
c
      do 20 i=1,nproc
        if(infoa(1,i).ne.iposn .or. infoa(2,i).ne.itt)goto 40
  20  continue
c
c  send proceed message to all processes
c
      call m_send(MSG_CONT)
      return
c
c  print information arrays every few seconds
c  use sleep routine to sleep (really we want to wait for a few
c  seconds before reprinting status arrays or until data arrives).
c
  40  continue
      call prininfo(4)
#ifndef cray-t3d
      call sleep(1)
#endif
      goto 10
      end

      subroutine mchkarch(iposn0)
c
c=======================================================================
c  master subroutine to checkpoint until if selected slaves have
c  reached position iposn on current timestep 
c  This routine is used by the master to verify that restart data
c  sent to selected slaves has been received successfully before
c  proceeding to process the next part of the restart dataset.
c
c  input:  iposn - position for checking
c=======================================================================
c
#include "param.h"
#include "ctmngr.h"
#include "mesdta.h"
c
c  set position
c
      iposn = iposn0
      call prininfo(3)
c
c  receive and process all waiting messages
c
  10  call m_recv(MSG_ANY)
c
c  check position and timestep correct for all processors
c
      do 20 ii=1,nsarch
        i=narcha(ii)
        if(infoa(1,i).ne.iposn .or. infoa(2,i).ne.itt)goto 40
  20  continue
c
c reset status position record ready for the next part of the restart
c data file. This ensures that no slave is sent the next part of
c the restart data before all slaves have indicated receipt of the
c previous part.
c
      do 25 ii=1,nsarch
       i=narcha(ii)
       infoa(1,i) = iposn + 1
  25  continue
c
      return
c
c  print information arrays every few seconds
c  use sleep routine to sleep (really we want to wait for a few
c  seconds before reprinting status arrays or until data arrives).
c
  40  continue
      call prininfo(4)
#ifndef cray-t3d
      call sleep(0)
#endif
      goto 10
      end
 
      subroutine inittsi
c
c=======================================================================
c  subroutine to initialise tso processing
c=======================================================================
c
#include "param.h"
#include "cdiag.h"
#include "ctmngr.h"
#include "mesdta.h"
#include "scalar.h"
c
      ltsip   = .true.
      itttsi  = itt
      timetsi = stamp
      volume = 0.
      ektot = 0.
      dtabs(1) = 0.
      dtabs(2) = 0.
      nstsi   = nproc
      do 10 i=1, nproc
        ltsia(i) = .false.
  10  continue
      call m_send(MSG_RQ_TSI)
      return
      end

      subroutine initsnap
c
c=======================================================================
c  subroutine to initialise snapshot processing
c=======================================================================
c
#include "param.h"
#include "ctmngr.h"
#include "mesdta.h"
c
      lsnapp = .true.
      ittsnap = itt
c
c Prepare to collect first variable
c snapreq will increment idsnp, so set to zero initially
c
      idsnp = 0
c
c------------------------------------------------------------------
c  open snapshot file name (call snapwr with idsnp = 0)
c------------------------------------------------------------------
c
      call snapwr
c
      return
      end

      subroutine initarch
c
c=======================================================================
c  subroutine to initialise archive processing
c=======================================================================
c
#include "param.h"
#include "ctmngr.h"
#include "mesdta.h"
c
      larchp = .true.
      larchq = .false.
      ittarch = itt
      do 20 j=1, NARCHV
       lidvra(j) = .false.
  20  continue
c
c Make copies of the time variables
c
      ittsav = itt
      totssav= totsec
      totdsav= totday
      yrssav = years
c
c Prepare to collect first variable
c archreq will increment idvar, so set to zero initially
c
      idvar = 0
      jarchu = 0
c
c------------------------------------------------------------------
c  create archive file name (call archwr with idvar = 0)
c------------------------------------------------------------------
c
      call archwr
c
      return
      end


      subroutine prininfo(iopt)
c
c=======================================================================
c  subroutine to print position information
c  if  iopt = 1 - reset main loop time flag (print after waiting statms s)
c           = 2 - check flags and if necessary print
c           = 3 - reset waiting loop time flag (print after waiting statis
c                 seconds with no messages arriving)
c           = 4 - as 2 but master at checkpoint
c=======================================================================
c
#include "param.h"
#include "mesdta.h"
#include "switch.h"
c
      integer tnow, time, tmain, twait, twaiti
      data tmain,twait /0,0/
      save tmain,twait,twaiti
c
c  find the current time
c
#ifdef cray-t3d
      tnow = rtc()*6.6e-9-time_start
#else
      tnow = time()-intime
#endif
c
c  entry points to reset clocks
c
      if(iopt.eq.1)then
        if(statms.gt.0.0)then
          tmain = (int(tnow/statms)+1)*statms
          tmain = max(tmain,tnow+1)
        else
          tmain = 0
        endif
        return
      elseif(iopt.eq.3)then
        twait = tnow
        if(statis.gt.0.0)then
          twaiti = max(1,int(statis))
        else
          twaiti= 0
        endif
        return
      endif
c
c  if the main loop status message is due or if no messages
c  have been received from the slaves for some time then
c  print the status arrays
c 
      if((statms.gt.0.0.and.tnow.gt.tmain).or.
     &   (twaiti.gt.0.and.tnow.ge.twait+twaiti))then
        now = mod(tnow,100000)
        if(iopt.eq.2)then
          write(stdout,31) now
        else
          write(stdout,32) now
        endif
        if(statms.gt.0.0.and.tnow.gt.tmain)then
          tmain = (int(tnow/statms)+1)*statms
          tmain = max(tmain,tnow+1)
        endif
        if(twaiti.gt.0.and.tnow.ge.twait+twaiti)then
          write(stdout,36)tnow-twait
          twaiti = twaiti+twaiti
        endif
        print 33
        do 40 i=1,nproc
        ii = infoa(1,i)
        if(ii.ge.0.and.ii.le.40)then
          mposmg = mesind(infoa(1,i))
        else
          mposmg = 0
        endif
        write(stdout,34)i,infoa(1,i),infoa(2,i),infoa(3,i),
     &                  messtr(mposmg)
  40    continue
        write(stdout,35)
        call flush(stdout)
      endif
  31  format(' ',/,' Status arrays at time =', i6,' s.',
     &             '  Master in main loop')
  32  format(' ',/,' Status arrays at time =', i6,' s.'
     &             '  Master waiting at checkpoint')
  33  format(' slave,  posn,  slave itt,  slave itbt,',
     &       '  Last reported: ')
  34  format(i5,i8,i10,i12,7x,a30)
  35  format(' ',/)
  36  format(' NOTE: Master has not received any messages from the',
     &       ' slaves for ',i5,' s.')
      return
      end


      subroutine m_send(msgtype)
c
c=======================================================================
c generic message handling subroutine for the OCCAM model.
c This subroutine is used by the master program to pack
c and send messages to the slave.
c=======================================================================
c
#include "def_master.h"
#include "param.h"
#include "scalar.h"
#include "switch.h"
#include "iounit.h"
#include "ctmngr.h"
#include "chmix.h"
#include "cvmix.h"
#include "coord.h"
#include "grdvar.h"
#include "levind.h"
#include "archive.h"
#include "snaps.h"
#include "varinfo.h"
#include "mesdta.h"
c
c=======================================================================
c pack and send command to abort (message MSG_ABORT)
c=======================================================================
c
      if(msgtype.eq.MSG_ABORT) then
        call pvmfinitsend( PVMDEFAULT, info )
        call pvmfmcast( nproc, itida, MSG_ABORT, info )
        if(idebug.gt.2)then
         print *,' Master:  Message MSG_ABORT sent to all slaves.'
        endif
        call mchkpnt(1)
        call closedown(info)
        stop
c
c=======================================================================
c pack and send command to continue (message MSG_CONT)
c=======================================================================
c
      elseif(msgtype.eq.MSG_CONT) then
        call pvmfinitsend( PVMDEFAULT, info )
        call pvmfmcast( nproc, itida, MSG_CONT, info )
        if(idebug.gt.2)then
         print *,' Master:  Message MSG_CONT sent to all slaves.'
        endif
c
c=======================================================================
c pack and send tid information from the master
c=======================================================================
c
      elseif(msgtype.eq.MSG_IDS) then
        call pvmfinitsend( PVMDEFAULT, info )
        call pvmfpack( INTEGER_PVM, nproc, 1, 1, info )
        call pvmfpack( INTEGER_PVM, itida, nproc, 1, info )
#ifdef cray-t3d
        call pvmfmcast( n$pes-1, itida, MSG_IDS, info )
#else
        call pvmfmcast( nproc, itida, MSG_IDS, info )
#endif
        if(idebug.gt.2)then
          write(stderr,*)' Master:  Message MSG_IDS sent to all slaves.'
        endif
c
c=======================================================================
c pack and send control information from the master
c=======================================================================
c
      elseif(msgtype.eq.MSG_CONTROL) then
        call pvmfinitsend( PVMDEFAULT, info )
        call pvmfpack( INTEGER_PVM,  init, 1, 1, info )
        call pvmfpack( INTEGER_PVM,restrt, 1, 1, info )
        call pvmfpack( INTEGER_PVM,    eb, 1, 1, info )
        call pvmfpack( REAL_PVM,     days, 1, 1, info )
        call pvmfpack( REAL_PVM,      tsi, 1, 1, info )
        call pvmfpack( REAL_PVM,     acor, 1, 1, info )
        call pvmfpack( REAL_PVM,   dgnstc, 1, 1, info )
        call pvmfpack( REAL_PVM,    snapd, 1, 1, info )
        call pvmfpack( REAL_PVM,    archd, 1, 1, info )
        call pvmfpack( REAL_PVM,     dtts, 1, 1, info )
        call pvmfpack( REAL_PVM,     dtuv, 1, 1, info )
        call pvmfpack( REAL_PVM,     dtbt, 1, 1, info )
        call pvmfpack( REAL_PVM,       am, 1, 1, info )
        call pvmfpack( REAL_PVM,       ah, 1, 1, info )
        call pvmfpack( REAL_PVM,     fkpm, 1, 1, info )
        call pvmfpack( REAL_PVM,     fkph, 1, 1, info )
        call pvmfpack( REAL_PVM,    cdbot, 1, 1, info )
        call pvmfpack( REAL_PVM,    dxdeg,    1, 1, info )
        call pvmfpack( REAL_PVM,    dydeg,    1, 1, info )
        call pvmfpack( REAL_PVM,       dz,   KM, 1, info )
#ifdef de_checkbd
        call pvmfpack( REAL_PVM,   dchkbd, 1, 1, info )
#endif
        call pvmfpack( INTEGER_PVM,  nmix, 1, 1, info )
        call pvmfpack( INTEGER_PVM,  ncon, 1, 1, info )
        call pvmfpack( INTEGER_PVM,  ntbt, 1, 1, info )
        call pvmfpack( INTEGER_PVM,idebug, 1, 1, info )
        call pvmfpack( INTEGER_PVM,iorest, 1, 1, info )
        call pvmfpack( INTEGER_PVM,ispvar, NSNAPS, 1, info )
        call pvmfpack( INTEGER_PVM,isplev, NSNAPS, 1, info )
        call pvmfpack( STRING, fnrest,80, 1, info )
c
        call pvmfmcast( nproc, itida, MSG_CONTROL, info )
        if(idebug.gt.2)then
          print *,' Master:  Message MSG_CONTROL sent to all slaves.'
        endif
c
c=======================================================================
c pack and send kmt and grid information
c=======================================================================
c
      elseif(msgtype.eq.MSG_TOPOG) then
       do 20 n=1,nproc
         if(idebug.gt.5)then
           print *,' Master:  Sending MSG_TOPOG to slave:',n
         endif
         call pvmfinitsend( PVMDEFAULT, info )
#ifdef pvm_buffer
         bufout(2) = it_l(n)
         bufout(3) = jt_l(n)
         bufout(4) = icycl_M(n)
         bufout(5) = iswm1_M(n)
         bufout(6) = jswm1_M(n)
         bufout(7) = swlo(n)
         bufout(8) = swla(n)
         lbuff = 8
#else
         call pvmfpack( INTEGER_PVM,  it_l(n),        1, 1, info )
         call pvmfpack( INTEGER_PVM,  jt_l(n),        1, 1, info )
         call pvmfpack( INTEGER_PVM,  icycl_M(n),        1, 1, info )
         call pvmfpack( INTEGER_PVM,  iswm1_M(n),        1, 1, info )
         call pvmfpack( INTEGER_PVM,  jswm1_M(n),        1, 1, info )
         call pvmfpack( REAL_PVM,     swlo(n),        1, 1, info )
         call pvmfpack( REAL_PVM,     swla(n),        1, 1, info )
#endif
         do 10 j=jt_s(n),jt_e(n)
         do 10 i=it_s(n),it_e(n)
          ii = mod(i-1,IMT_M) + 1
#ifdef pvm_buffer
             lbuff = lbuff+1
             if(lbuff+1.gt.LENBUF)then
               print *, ' error T1 in m_send while packing',
     &             ' MSG_TOPOG for slave: ',n,itida(n)
               print *, ' buffer overflow: lbuff =',lbuff,
     &                    ' ii = ',ii,' j = ',j
             endif
             bufout(lbuff) = kmt(ii,j)
#else
           call pvmfpack(INTEGER_PVM , kmt(ii,j),  1, 1, info )
#endif
 10      continue
         do 15 j=jt_s(n),jt_e(n)
         do 15 i=it_s(n),it_e(n)
           ii = mod(i-1,IMT_M) + 1
#ifdef pvm_buffer
             lbuff = lbuff+1
             if(lbuff+1.gt.LENBUF)then
               print *, ' error T2 in m_send while packing',
     &             ' MSG_TOPOG for slave: ',n,itida(n)
               print *, ' buffer overflow:  lbuff =',lbuff,
     &                    ' ii = ',ii,' j = ',j
             endif
             bufout(lbuff) = kmu(ii,j)
#else
           call pvmfpack(INTEGER_PVM , kmu(ii,j),  1, 1, info )
#endif
 15      continue
         do 18 j=jt_s(n),jt_e(n)
         do 18 i=it_s(n),it_e(n)
           ii = mod(i-1,IMT_M) + 1
#ifdef pvm_buffer
             lbuff = lbuff+1
             if(lbuff+1.gt.LENBUF)then
               print *, ' error T3 in m_send while packing',
     &             ' MSG_TOPOG for slave: ',n,itida(n)
               print *, ' buffer overflow: lbuff =',lbuff,
     &                    ' ii = ',ii,' j = ',j
             endif
             bufout(lbuff) = kpn(ii,j)
#else
           call pvmfpack(INTEGER_PVM , kpn(ii,j),  1, 1, info )
#endif
 18        continue
#ifdef pvm_buffer
         bufout(1) = lbuff
         call pvmfpack(REAL_PVM, bufout, lbuff, 1, info)
#endif
c
         call pvmfsend( itida(n), msgtype, info )
 20      continue
       if(idebug.gt.2)then
         print *,' Master:  Message MSG_TOPOG sent to all slaves.'
       endif
c
c=======================================================================
c  entry point to send restart time variables to all slaves
c=======================================================================
c
      elseif(msgtype.eq.MSG_TIMEVAR) then
        call pvmfinitsend( PVMDEFAULT, info )
        call pvmfpack( INTEGER_PVM, itt   , 1, 1, info )
        call pvmfpack( REAL_PVM  , totsec, 1, 1, info )
        call pvmfpack( REAL_PVM  , totday, 1, 1, info )
        call pvmfpack( REAL_PVM  , years , 1, 1, info )
        call pvmfmcast(nproc, itida, MSG_TIMEVAR, info )
        if(idebug.gt.2)then
          print *,
     &    ' Master:  Message MSG_TIMEVAR sent to all slaves.'
        endif
c
c=======================================================================
c  entry point to send restart data to the slaves 
c=======================================================================
c
      elseif(msgtype.eq.MSG_RESTART) then
        do 50 nn=1,nsarch
         n=narcha(nn)
         jlow=max(jarchl,jt_s(n))
         jupp=min(jarchu,jt_e(n))
         call pvmfinitsend( PVMDEFAULT, info )
#ifdef pvm_buffer
           bufout(2) = idvar
           bufout(3) = jlow
           bufout(4) = jupp
           lbuff = 4
#else
         call pvmfpack( INTEGER_PVM, idvar, 1, 1, info )
         call pvmfpack( INTEGER_PVM, jlow , 1, 1, info )
         call pvmfpack( INTEGER_PVM, jupp , 1, 1, info )
#endif
c
         if(idvar.le.NUM2D) then
          do 60 j=jlow,jupp
          do 60 i=it_s(n),it_e(n)
           ii = mod(i-1,IMT_M) + 1
#ifdef pvm_buffer
           lbuff = lbuff+1
           if(lbuff+1.gt.LENBUF)then
               print *, ' error in m_send while packing',
     &             ' MSG_RESTART for slave: ',n,itida(n)
               print *, ' buffer overflow:  lbuff =',lbuff
           endif
           bufout(lbuff) = rest2d(ii,j)
#else
           call pvmfpack(REAL_PVM, rest2d(ii,j), 1, 1, info)
#endif
60        continue
c
         elseif(vartgrd(idvar)) then
          joff=((jarchu-1)/JSUB_M)*JSUB_M
          do 70 j=jlow,jupp 
           jj=j-joff
          do 70 i=it_s(n),it_e(n)
           ii = mod(i-1,IMT_M) + 1
           if(kmt(ii,j).ne.0) then
#ifdef pvm_buffer
             if(lbuff+kmt(ii,j)+1.gt.LENBUF)then
               print *, ' error in m_send while packing',
     &             ' MSG_RESTART for slave: ',n,itida(n)
               print *, ' buffer overflow:  lbuff =',lbuff
             endif
             do 61 k=1,kmt(ii,j)
               lbuff = lbuff+1
               bufout(lbuff) = rest3d(k,ii,jj)
  61         continue
#else
             call pvmfpack(REAL_PVM,  
     &                        rest3d(1,ii,jj),kmt(ii,j), 1, info)
#endif
           endif 
70        continue
         else
c
          joff=((jarchu-1)/JSUB_M)*JSUB_M
          do 80 j=jlow,jupp 
           jj=j-joff
          do 80 i=it_s(n),it_e(n)
           ii = mod(i-1,IMT_M) + 1
           if(kmu(ii,j).ne.0) then
#ifdef pvm_buffer
             if(lbuff+kmu(ii,j)+1.gt.LENBUF)then
               print *, ' error in m_send while packing',
     &             ' MSG_RESTART for slave: ',n,itida(n)
               print *, ' buffer overflow:  lbuff =',lbuff
             endif
             do 62 k=1,kmu(ii,j)
               lbuff = lbuff+1
               bufout(lbuff) = rest3d(k,ii,jj)
  62         continue
#else
             call pvmfpack(REAL_PVM,  
     &                        rest3d(1,ii,jj),kmu(ii,j), 1, info)
#endif
            endif 
80        continue
         endif
#ifdef pvm_buffer
         bufout(1) = lbuff
         call pvmfpack(REAL_PVM, bufout, lbuff, 1, info)
#endif
         call pvmfsend(iarcha(nn), MSG_RESTART, info )
50      continue
        if(idebug.gt.2)then
          print *,' Master:  Message MSG_RESTART sent to',
     &           ' selected slaves. idvar = ',idvar
        endif
c
c=======================================================================
c  entry point to ask for tsi information (message MSG_RQ_TSI)
c  from all slaves for timestep itttsi 
c=======================================================================
c
      elseif(msgtype.eq.MSG_RQ_TSI) then
        call pvmfinitsend( PVMDEFAULT, info )
        call pvmfpack( INTEGER_PVM, itttsi, 1, 1, info )
        call pvmfmcast(nproc, itida, MSG_RQ_TSI, info )
        if(idebug.gt.2)then
          print *,' Master:  Message MSG_RQ_TSI sent to all slaves.'
        endif
c
c=======================================================================
c  entry point to ask for part of the archive data (message MSG_RQ_ARC)
c  from selected slaves for timestep ittarch
c=======================================================================
c
      elseif(msgtype.eq.MSG_RQ_ARC) then
        call pvmfinitsend( PVMDEFAULT, info )
        call pvmfpack( INTEGER_PVM, ittarch, 1, 1, info )
        call pvmfpack( INTEGER_PVM, idvar  , 1, 1, info )
        call pvmfpack( INTEGER_PVM, jarchl , 1, 1, info )
        call pvmfpack( INTEGER_PVM, jarchu , 1, 1, info )
        call pvmfmcast(nsarch, iarcha, MSG_RQ_ARC, info )
        if(idebug.gt.2)then
          print *,
     &    ' Master:  Message MSG_RQ_ARC sent to selected slaves.',
     &    ' idvar, range = ',idvar,jarchl,jarchu
        endif
c
c=======================================================================
c  entry point to inform the slaves that all archive data 
c  has been received and that they may clear their buffers (MSG_ARC_CLR)
c=======================================================================
c 
      elseif(msgtype.eq.MSG_ARC_CLR) then 
        call pvmfinitsend( PVMDEFAULT, info ) 
        call pvmfpack( INTEGER_PVM, ittarch, 1, 1, info )   
        call pvmfmcast(nproc, itida, MSG_ARC_CLR, info ) 
        if(idebug.gt.2)then 
          print *, 
     &    ' Master:  Message MSG_ARC_CLR sent to all slaves.' 
        endif 
c
c=======================================================================
c  entry point to send snapshot request to the slaves
c=======================================================================
c
      elseif(msgtype.eq.MSG_RQ_SNAP) then
        call pvmfinitsend( PVMDEFAULT, info )
        call pvmfpack( INTEGER_PVM,ittsnap, 1, 1, info )
        call pvmfpack( INTEGER_PVM,  idsnp, 1, 1, info )
        call pvmfmcast(nproc, itida, MSG_RQ_SNAP, info )
        if(idebug.gt.2)then
          print *,
     &    ' Master:  Message MSG_RQ_SNAP sent to all slaves.'
        endif
c
c=======================================================================
c  unrecognised request
c=======================================================================
c
      else
        write(stderr,*) 
     &  ' Master:  Message type: ',msgtype,' not yet implemented'
      endif
c
c=======================================================================
c  flush stderr to ensure debug messages are received
c=======================================================================
c
      if(idebug.gt.2)then
          call flush(stderr)
      endif
      return
      end


      subroutine m_recv(msgtype)
c
c=======================================================================
c   Generic message handling subroutine for the master.
c   This subroutine is used by the master program to receive,
c   unpack and, if necessary, process messages.
c
c   if msgtype .ge. 0 receive only requested message then return
c                     otherwise receive messages until none left
c
c  NOTE:  I do not think the >0 option is essential.  If it is not
c         used, delete it later.  (djw)
c=======================================================================
c
#include "def_master.h"
#include "param.h"
#include "scalar.h"
#include "switch.h"
#include "iounit.h"
#include "cdiag.h"
#include "chmix.h"
#include "levind.h"
#include "archive.h"
#include "snaps.h"
#include "varinfo.h"
#include "mesdta.h"
c
      integer msgtype, time
c
c=======================================================================
c  check to see if anything to receive
c  return if no messages
c=======================================================================
c
  10  if(msgtype.ge.0)then
        call pvmfrecv (-1, msgtype, ibuff)
      else
        call pvmfnrecv (-1, MSG_ANY, ibuff)
        if(ibuff.eq.0)return
      endif
c
c  message received - reset waiting flag
c
      call prininfo(3)
c
c=======================================================================
c  get information about the buffer
c=======================================================================
c
      call pvmfbufinfo(ibuff,ibytes, msgtag, itid, info)
#ifdef cray-t3d
        call pvmfgetpe(itid,itidpe)
        itid = itidpe
#endif
c
c  read name of slave and check it corresponds to the tid
c
      call pvmfunpack( INTEGER_PVM, mes,  1, 1, info )
      if(mes.lt.1.or.mes.gt.MXSLAVE)then
        write(stderr,*) ' Error: Message ',msgtag,' from tid ',itid,
     &     ' has slave number ',mes,' out of range 0 - ',MXSLAVE
      elseif(itid.ne.itida(mes))then
        write(6,*) ' Error: Message ',msgtag,' from slave ',mes,
     &     ' has tid ',itid,', should be ',itida(mes)
      endif
c
c=======================================================================
c  process abort request message (MSG_BYE)
c=======================================================================
c
      if(msgtag.eq.MSG_BYE)then
        write(stderr,*) ' Abort request received from slave: ',mes
        write(stderr,*) ' Closing down all slaves'
        call m_send(MSG_ABORT)
c
c=======================================================================
c  process position message (MSG_POSN)
c  message contains:  me, itt, itt0, iposn
c=======================================================================
c
      elseif(msgtag.eq.MSG_POSN)then
        call pvmfunpack( INTEGER_PVM, infoa(1,mes),  1, 1, info )
        call pvmfunpack( INTEGER_PVM, infoa(2,mes),  1, 1, info )
        call pvmfunpack( INTEGER_PVM, infoa(3,mes),  1, 1, info )
        if(idebug.gt.2)then
          write(stderr,*) ' Master:  Message MSG_POSN received: ',
     &    mes, infoa(1,mes),infoa(2,mes),infoa(3,mes)
        endif
c
c=======================================================================
c process information message (MSG_OUTSTR)
c message contains: me, string for printing
c=======================================================================
c
      elseif(msgtag.eq.MSG_OUTSTR) then
        call pvmfunpack( STRING,   outstr,80, 1, info )
        if(idebug.gt.2)then
          write(stderr,*) ' Master:  Message MSG_OUTSTR received ',
     &                      'from slave ',mes,' ',outstr
        else
          write(stderr,921)mes,outstr
 921  format(' msg from',i4,':  ',a)
        endif
c
c=======================================================================
c  process timestep 'tsi' information (MSG_RY_TSI)
c=======================================================================
c
      elseif(msgtag.eq.MSG_RY_TSI)then
        call pvmfunpack( INTEGER_PVM,itts,    1, 1, info )
        call pvmfunpack( REAL_PVM,   vols,    1, 1, info )
        call pvmfunpack( REAL_PVM, ektots,    1, 1, info )
        call pvmfunpack( REAL_PVM, dtemps,    1, 1, info )
        call pvmfunpack( REAL_PVM, dsalts,    1, 1, info )
c
c  checks
c
        if(.not.ltsip)
     &    write(stderr,*) ' Error: Message MSG_RY_TSI from slave ',mes,
     &                  ' sent when tsi processing inactive.' 
        if(itts.ne.itttsi)then
          write(stderr,*) ' Error: Message MSG_RY_TSI from slave ',mes,
     &                    ' has wrong time flag'
          write(stderr,*) '     : itt from slave is ',itts,
     &                    ' should be ', itttsi
        endif          
        if(ltsia(mes))
     &    write(stderr,*) ' Error: Message MSG_RY_TSI from slave ',mes,
     &                    ' sent twice'
        if(idebug.gt.2)then
          write(stderr,*) ' Master:  Message MSG_RY_TSI received: ',
     &                      mes,vols,ektots,dtemps,dsalts
        endif
c
        ltsia(mes) = .true.
        nstsi = nstsi - 1
        volume   = volume + vols
        ektot    = ektot + ektots
        dtabs(1) = dtabs(1) + dtemps
        dtabs(2) = dtabs(2) + dsalts
        if(idebug.gt.2)then
        print *,' Master:  Message MSG_RY_TSI. volume, ektot,dtabs =',
     &            volume, ektot, dtabs(1), dtabs(2)
        endif
c
c print tsi line
c
        if(nstsi.eq.0.and.ltsip)then
          ektot    = ektot/volume
          dtabs(1) = dtabs(1)/volume
          dtabs(2) = dtabs(2)/volume
c
#ifdef cray-t3d
          total_time=rtc()*6.6e-9-time_start
          write(stderr,9602) itts,timetsi,ektot,dtabs(1),dtabs(2),
     &                       total_time
9602      format(1x,'ts=',i7,1x,a32,' ke=', 1pe13.6,
     &    ' dtemp=',1pe13.6,' dsalt=',1pe13.6,' time=',0pf10.4)
#else
          itime = time()-intime
          write(stderr,9602) itts,timetsi,ektot,dtabs(1),dtabs(2),
     &                       itime
9602      format(1x,'ts=',i7,1x,a32,' ke=', 1pe13.6,
     &    ' dtemp=',1pe13.6,' dsalt=',1pe13.6,' time=',i6)
#endif
          call flush(stderr)
          ltsip = .false.
        endif
c
c=======================================================================
c  process archive replies (MSG_RY_ARC)
c=======================================================================
c
      elseif(msgtag.eq.MSG_RY_ARC)then
#ifdef pvm_buffer
        call pvmfunpack (REAL_PVM, bufout(1), 1, 1, info)
        lbuff = nint(bufout(1))
        if(lbuff.gt.LENBUF)then
           write(stderr,*) ' error in m_recv when receiving ',
     &                      ' MSG_RY_SNAP'
           write(stderr,*) ' buffer overflow'
        endif
        call pvmfunpack (REAL_PVM, bufout(2), lbuff-1, 1, info)
        ittsav  = nint(bufout(2))
        idvar2  = nint(bufout(3))
        lbuff   = 3
#else
        call pvmfunpack(INTEGER_PVM,ittsav,    1, 1, info )
        call pvmfunpack(INTEGER_PVM,idvar2, 1, 1, info)
#endif
c
c check for expected timestep
c
        if(.not.larchp) write(stderr,*)
     &                  ' Error: Message MSG_RY_ARC from slave ',mes,
     &                  ' sent when archive processing inactive.'
        if(ittsav.ne.ittarch) then
          write(stderr,*) ' Error: Message MSG_RY_ARC from slave ',mes,
     &                    ' has wrong time flag'
          write(stderr,*) '     : itt from slave is ',ittsav,
     &                    ' should be ', ittarch
        endif
c
c
c check this is the expected variable type
c
       if(.not.lidvra(idvar2))write(stderr,*)     
     &                ' Error: Message MSG_RY_ARC from slave ',mes,
     &                ' contains wrong variable type'
#ifdef pvm_buffer
 55    imast = nint(bufout(lbuff+1))
#else
 55    call pvmfunpack(INTEGER_PVM,imast, 1, 1, info)
#endif
c
c check for "end of message" marker
c
       if(imast.gt.0) then
#ifdef pvm_buffer
         lbuff = lbuff+2
         jmast = nint(bufout(lbuff))
#else
         call pvmfunpack(INTEGER_PVM, jmast, 1, 1, info)
#endif
c
c check for corrupted message
c
          if(imast.gt.IMT_M.or.jmast.gt.JMT_M) then
            write(stderr,*) ' Invalid master indices: ',imast,jmast
            write(stderr,*) ' currently unpacking ',nseapt,
     &                      'th sea-point'
            call m_send(MSG_ABORT)
          endif
c
c check this reply is unique
c
          if(lregra(imast,jmast)) then
            write(stderr,*) ' Error: Message MSG_RY_ARC from slave ',
     &                   mes, ' is overwriting previously received ',
     &                        ' data or land points'
          else
c
c Increment a sea-point counter
c
            nseapt = nseapt + 1
            lregra(imast,jmast) = .true.
          endif
c 
c Now unpack requested variable(s) for this point 
c 
          if(idvar2.le.NUM2D) then 
#ifdef pvm_buffer
            lbuff = lbuff+1
            rest2d(imast,jmast) = bufout(lbuff)
#else
            call pvmfunpack(REAL_PVM, rest2d(imast,jmast), 1, 1, info) 
#endif
          elseif(vartgrd(idvar2)) then
c
c convert jmast to sub-volume index
c
            jzon = jmast - ((jmast-1)/JSUB_M)*JSUB_M
#ifdef pvm_buffer
            do 56 k=1,kmt(imast,jmast)
              lbuff = lbuff+1
              rest3d(k,imast,jzon) = bufout(lbuff)
  56        continue
#else
            call pvmfunpack(REAL_PVM, rest3d(1,imast,jzon), 
     &                               kmt(imast,jmast), 1, info)
#endif
          else
c
c convert jmast to sub-volume index
c
           jzon = jmast - ((jmast-1)/JSUB_M)*JSUB_M 
#ifdef pvm_buffer
           do 57 k=1,kmu(imast,jmast)
             lbuff = lbuff+1
             rest3d(k,imast,jzon) = bufout(lbuff)
  57       continue
#else
           call pvmfunpack(REAL_PVM, rest3d(1,imast,jzon), 
     &                               kmu(imast,jmast), 1, info)
#endif
         endif   
         goto 55
       endif
c
       if(nseapt.ge.nseachk) call archwr       
c
       if(idebug.gt.2)then
          print *,' Master:  Message MSG_RY_ARC from slave ',
     &            mes,'  idvar, itt, nseapt, nseachk = ',
     &                   idvar,ittsav,nseapt,nseachk
       endif
c
c=======================================================================
c  process snapshot replies (MSG_RY_SNAP)
c=======================================================================
c
      elseif(msgtag.eq.MSG_RY_SNAP)then
#ifdef pvm_buffer
        call pvmfunpack (REAL_PVM, bufout(1), 1, 1, info)
        lbuff = nint(bufout(1))
        if(lbuff.gt.LENBUF)then
           write(stderr,*) ' error in m_recv when receiving ',
     &                      ' MSG_RY_SNAP'
           write(stderr,*) ' buffer overflow'
        endif
        call pvmfunpack (REAL_PVM, bufout(2), lbuff-1, 1, info)
        ittsnap2  = nint(bufout(2))
        idsnp2    = nint(bufout(3))
        lbuff     = 3
#else
        call pvmfunpack(INTEGER_PVM,ittsnap2,    1, 1, info )
        call pvmfunpack(INTEGER_PVM, idsnp2, 1, 1, info)
#endif
c
c check for expected timestep and variable type
c
        if(.not.lsnapp) then 
          write(stderr,*)' Error: Message MSG_RY_SNAP from slave '
     &              ,mes,' sent when snapshot processing inactive.'
        endif
        if(ittsnap2.ne.ittsnap) then 
          write(stderr,*) ' Error: Message MSG_RY_SNAP from slave ',
     &                      mes,' has wrong time flag'
          write(stderr,*) '      : itt from slave is ',ittsnap2,
     &                    ' should be ', ittsnap
        endif
        if(idsnp2.ne.idsnp) then
          write(stderr,*) ' Error: Message MSG_RY_SNAP from slave ',
     &                      mes,' contains wrong variable type'
        endif
#ifdef pvm_buffer
 65      imast = nint(bufout(lbuff+1))
#else
 65      call pvmfunpack(INTEGER_PVM, imast, 1, 1, info)
#endif
         if(imast.gt.0)then
#ifdef pvm_buffer
           jmast = nint(bufout(lbuff+2))
#else
           call pvmfunpack(INTEGER_PVM, jmast, 1, 1, info)
#endif
c
c check for corrupted message
c
          if(imast.gt.IMT_M.or.jmast.gt.JMT_M) then
            write(stderr,*) ' Invalid master indices: ',imast,jmast
            write(stderr,*) ' currently unpacking ',idsnp,'th snapshot'
            call m_send(MSG_ABORT)
          endif
#ifdef pvm_buffer
          snap2d(imast,jmast) = bufout(lbuff+3)
          lbuff = lbuff+3
#else
          call pvmfunpack(REAL_PVM, snap2d(imast,jmast), 1, 1, info)
#endif
          goto 65
        endif
c
c Decrement the reply counter
c
       nssnap = nssnap - 1
c
       if(nssnap.eq.0) call snapwr
c
        if(idebug.gt.2)then
          print *,' Master:  Message MSG_RY_SNAP from slave ',
     &            mes,'  idsnp, itt, nssnap = ',
     &                   idsnp2,ittsnap2,nssnap         
        endif
c
c=======================================================================
c  error          
c=======================================================================
c
      else
        write(stderr,*) ' Master:  Message type: ',msgtag,
     &                  ' not yet implemented'
      endif
c
c=======================================================================
c  flush stderr to ensure debug messages are received
c=======================================================================
c
      if(idebug.gt.2)then
          call flush(stderr)
      endif
      if(msgtype.ne.msgtag)goto 10
      return
      end
