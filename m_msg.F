      subroutine m_initmess
c
c=======================================================================
c @(#) SCCS module: m_msg.F  version: 1.1
c     Creation date: 10/13/97
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
        lmeta(i)  = .true.
  20  continue
      lsend  = .true.
      ltsip  = .false.
      lsnapp = .false.
      lsnapq = .false.
      larchp = .false.
      larchq = .false.
      lmetp  = .true.
      lmetq  = .false.
c
      return
      end

      subroutine m_abort
c
c=======================================================================
c  subroutine to send abort message to slaves
c=======================================================================
c
#include "param.h"
#include "mesdta.h"
      call m_send(MSG_ABORT)
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
 
      subroutine nextmet
#ifdef hcomments
c
c @(#) SCCS module: m_msg.F  version: 1.1
c    Creation date: 10/13/97
c
c=======================================================================
c  subroutine to read in the next set of meteorological data and
c  send MSG_MET_RDY message
c=======================================================================
c
c 
#endif
#include "def_master.h"
#include "param.h"
#include "mesdta.h"
#include "ctmngr.h"
c
      call metrd(mnextp)
      nmeta = nproc
      lmetp = .true.
      call m_send(MSG_MET_RDY)
      return
      end

      subroutine inittsi
c
c=======================================================================
c  subroutine to initialise tsi processing
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
      tnow = rtc()*6.6e-9-start_time
#else
      tnow = MPI_WTIME()-start_time
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
#include "cvbc.h"
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
      MPII = MPI_INTEGER
      MPIR = MPI_REAL
      MPIL = LENBUF*4
      MPIC = MPI_COMM_WORLD
c
      if(msgtype.eq.MSG_ABORT) then
        do iproc=1,nproc
 100      call MPI_BSEND(bufout,1,MPII,iproc,msgtype,MPIC,ierr)
          call M_BERR(ierr,*100)
        enddo
        if(idebug.gt.2)then
          write(stderr,*)
     &         ' Master:  Message MSG_ABORT sent to all slaves.'
        endif
        lsend = .false.
        call mchkpnt(1)
        call closedown(info)
        stop
c
c=======================================================================
c  1.  Send message type MSG_CONT
c      command to continue
c=======================================================================
c
      elseif(msgtype.eq.MSG_CONT) then
        do iproc=1,nproc
 101      call MPI_BSEND(bufout,1,MPII,iproc,msgtype,MPIC,ierr)
          call M_BERR(ierr,*101)
        enddo
        if(idebug.gt.2)then
          write(stderr,'(a)')
     &         ' Master:  Message MSG_CONT sent to all slaves.'
        endif
c
c=======================================================================
c 2.  Send message type MSG_IDS
c     pack and send tid information to the slaves
c=======================================================================
c
      elseif(msgtype.eq.MSG_IDS) then
        do iproc=1,nproc
 102      call MPI_BSEND(nproc,1,MPII,iproc,msgtype,MPIC,ierr)
          call M_BERR(ierr,*102)
        enddo
        if(idebug.gt.2)then
          write(stderr,'(a)')
     &         ' Master:  Message MSG_IDS sent to all slaves.'
        endif
c
c=======================================================================
c  3.  Send message type MSG_CONTROL
c      pack and send control information to the slaves
c=======================================================================
c
      elseif(msgtype.eq.MSG_CONTROL) then
        ib = 0
        call MPI_PACK(init  ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(restrt,1,MPII,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(eb    ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(days  ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(tsi   ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(acor  ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(dgnsts,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(snapd ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(archd ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(dtts  ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(dtuv  ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(dtbt  ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(am    ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(ah    ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(fkpm  ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(fkph  ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(cdbot ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(dxdeg ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(dydeg ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(dz    ,KM,MPIR,bufout,MPIL,ib,MPIC,ierr)
#ifdef de_checkbd
        call MPI_PACK(dchkbd,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
#endif
        call MPI_PACK(nmix  ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(ncon  ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(ntbt  ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(idebug  ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(iorest  ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(ispvar  ,NSNAPS,MPII,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(isplev  ,NSNAPS,MPII,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(fnrest  ,80,MPI_CHARACTER,bufout,
     &                                            MPIL,ib,MPIC,ierr)
c
        do iproc=1,nproc
 103      call MPI_BSEND(bufout,ib,MPI_PACKED,iproc,msgtype,MPIC,ierr)
          call M_BERR(ierr,*103)
        end do
        if(idebug.gt.2)then
          write(stderr,'(a)')
     &         ' Master:  Message MSG_CONTROL sent to all slaves.'
        endif
c
c=======================================================================
c  4.  Send message type MSG_TOPOG
c      pack and send kmt and grid information
c=======================================================================
c
      elseif(msgtype.eq.MSG_TOPOG) then
       do 20 n=1,nproc
         if(idebug.gt.5)then
           print *,' Master:  Sending MSG_TOPOG to slave:',n
         endif
         bufout(2) = it_l(n)
         bufout(3) = jt_l(n)
         bufout(4) = icycl_M(n)
         bufout(5) = iswm1_M(n)
         bufout(6) = jswm1_M(n)
         bufout(7) = swlo(n)
         bufout(8) = swla(n)
         lbuff = 8
         do 10 j=jt_s(n),jt_e(n)
         do 10 i=it_s(n),it_e(n)
           ii = mod(i-1,IMT_M) + 1
           lbuff = lbuff+1
           if(lbuff+1.gt.LENBUF)then
             print *, ' error T1 in m_send while packing',
     &             ' MSG_TOPOG for slave: ',n
             print *, ' buffer overflow: lbuff =',lbuff,
     &                    ' ii = ',ii,' j = ',j
           endif
           bufout(lbuff) = kmt(ii,j)
 10      continue
         do 15 j=jt_s(n),jt_e(n)
         do 15 i=it_s(n),it_e(n)
           ii = mod(i-1,IMT_M) + 1
           lbuff = lbuff+1
           if(lbuff+1.gt.LENBUF)then
               print *, ' error T2 in m_send while packing',
     &             ' MSG_TOPOG for slave: ',n
               print *, ' buffer overflow:  lbuff =',lbuff,
     &                    ' ii = ',ii,' j = ',j
           endif
           bufout(lbuff) = kmu(ii,j)
 15      continue
         do 18 j=jt_s(n),jt_e(n)
         do 18 i=it_s(n),it_e(n)
           ii = mod(i-1,IMT_M) + 1
             lbuff = lbuff+1
             if(lbuff+1.gt.LENBUF)then
               print *, ' error T3 in m_send while packing',
     &             ' MSG_TOPOG for slave: ',n
               print *, ' buffer overflow: lbuff =',lbuff,
     &                    ' ii = ',ii,' j = ',j
             endif
             bufout(lbuff) = kpn(ii,j)
 18        continue
         bufout(1) = lbuff
c
 104     call MPI_BSEND(bufout,lbuff,MPI_REAL,
     &                          n,MSG_TOPOG,MPIC,ierr)
         call M_BERR(ierr,*104)
 20    continue
       if(idebug.gt.2)then
          write(stderr,'(a)')
     &    ' Master:  Message MSG_TOPOG sent to all slaves.'
       endif
c
c=======================================================================
c  entry point to send restart time variables to all slaves
c=======================================================================
c
      elseif(msgtype.eq.MSG_TIMEVAR) then
        bufout(1) = itt
        bufout(2) = totsec
        bufout(3) = totday
        bufout(4) = years
        do iproc=1,nproc
 105      call MPI_BSEND(bufout,4,MPI_REAL,
     &                          iproc,MSG_TIMEVAR,MPIC,ierr)
          call M_BERR(ierr,*105)
        enddo
        if(idebug.gt.2)then
          write(stderr,'(a)')
     &         ' Master:  Message MSG_TIMEVAR sent to all slaves.'
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
          bufout(2) = idvar
          bufout(3) = jlow
          bufout(4) = jupp
          lbuff = 4
c
          if(idvar.le.NUM2D) then
            do 60 j=jlow,jupp
            do 60 i=it_s(n),it_e(n)
              ii = mod(i-1,IMT_M) + 1
              lbuff = lbuff+1
              if(lbuff+1.gt.LENBUF)then
                print *, ' error in m_send while packing',
     &              ' MSG_RESTART for slave: ',n
                print *, ' buffer overflow:  lbuff =',lbuff
              endif
              bufout(lbuff) = rest2d(ii,j)
60          continue
c
          elseif(vartgrd(idvar)) then
            joff=((jarchu-1)/JSUB_M)*JSUB_M
            do 70 j=jlow,jupp 
            jj=j-joff
            do 70 i=it_s(n),it_e(n)
              ii = mod(i-1,IMT_M) + 1
              if(kmt(ii,j).ne.0) then
                if(lbuff+kmt(ii,j)+1.gt.LENBUF)then
                  print *, ' error in m_send while packing',
     &              ' MSG_RESTART for slave: ',n
                  print *, ' buffer overflow:  lbuff =',lbuff
                endif
                do 61 k=1,kmt(ii,j)
                  lbuff = lbuff+1
                  bufout(lbuff) = rest3d(k,ii,jj)
  61            continue
              endif 
70          continue
c
          else
            joff=((jarchu-1)/JSUB_M)*JSUB_M
            do 80 j=jlow,jupp 
            jj=j-joff
            do 80 i=it_s(n),it_e(n)
              ii = mod(i-1,IMT_M) + 1
              if(kmu(ii,j).ne.0) then
                if(lbuff+kmu(ii,j)+1.gt.LENBUF)then
                  print *, ' error in m_send while packing',
     &               ' MSG_RESTART for slave: ',n
                  print *, ' buffer overflow:  lbuff =',lbuff
                endif
                do 62 k=1,kmu(ii,j)
                  lbuff = lbuff+1
                  bufout(lbuff) = rest3d(k,ii,jj)
  62            continue
              endif 
80          continue
          endif
          bufout(1) = lbuff
c
c send data
c
 106      call MPI_BSEND(bufout,lbuff,MPI_REAL,iarcha(nn),
     &                          MSG_RESTART,MPIC,ierr)
          call M_BERR(ierr,*106)
50      continue
        if(idebug.gt.2)then
          write(stderr,'(a,a,i5)')
     &      ' Master:  Message MSG_RESTART sent to selected slaves.',
     &      ' idvar = ',idvar
        endif
c
c=======================================================================
c 10.  Send message type MSG_MET1 to slaves
c      initial met data including month one
c=======================================================================
c
        elseif(msgtype.eq.MSG_MET1)then
c  loop over slaves, clearing the input buffers on the way
          do 130 n=1,nproc
            call m_recv(MSG_ANY)
c
            ib = 0
            call MPI_PACK(mlast    ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
            call MPI_PACK(iw_l(n)  ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
            call MPI_PACK(jw_l(n)  ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
            call MPI_PACK(i_w      ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
            call MPI_PACK(swlo_w(n),1,MPIR,bufout,MPIL,ib,MPIC,ierr)
            call MPI_PACK(swla_w(n),1,MPIR,bufout,MPIL,ib,MPIC,ierr)
            call MPI_PACK(dx_w     ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
            call MPI_PACK(dy_w     ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
c
            do 120 j=jw_s(n),jw_e(n)
            do 120 i=iw_s(n),iw_e(n)
              ii=mod(i-1,NWX_M) + 1
              call MPI_PACK(wstx(ii,j),1,MPIR,bufout,MPIL,ib,MPIC,ierr)
              call MPI_PACK(wsty(ii,j),1,MPIR,bufout,MPIL,ib,MPIC,ierr)
 120        continue
c
            call MPI_PACK(is_l(n)  ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
            call MPI_PACK(js_l(n)  ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
            call MPI_PACK(is       ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
            call MPI_PACK(swlo_s(n),1,MPIR,bufout,MPIL,ib,MPIC,ierr)
            call MPI_PACK(swla_s(n),1,MPIR,bufout,MPIL,ib,MPIC,ierr)
            call MPI_PACK(dx_s     ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
            call MPI_PACK(dy_s     ,1,MPIR,bufout,MPIL,ib,MPIC,ierr)
c
            do 125 j=js_s(n),js_e(n)
            do 125 i=is_s(n),is_e(n)
              ii=mod(i-1,NSX_M) + 1
              call MPI_PACK(ssta(ii,j),1,MPIR,bufout,MPIL,ib,MPIC,ierr)
              call MPI_PACK(sssa(ii,j),1,MPIR,bufout,MPIL,ib,MPIC,ierr)
 125        continue
c
 107      call MPI_BSEND(bufout,ib,MPI_PACKED,n,
     &                          MSG_MET1,MPIC,ierr)
          call M_BERR(ierr,*107)
 130     continue
c
       if(idebug.gt.2)then
         write(stderr,'(a)')
     &        ' Master:  Message MSG_MET1 sent to all slaves.'
       endif
c
c=======================================================================
c 11.  Send message type MSG_MET2 to slaves
c=======================================================================
c
        elseif(msgtype.eq.MSG_MET2)then
c  loop over slaves, clearing the input buffers on the way
          do 160 n=1,nproc
            call m_recv(MSG_ANY)
c
            ib = 0
            call MPI_PACK(mnext    ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
            do 140 j=jw_s(n),jw_e(n)
            do 140 i=iw_s(n),iw_e(n)
              ii=mod(i-1,NWX_M) + 1
              call MPI_PACK(wstx(ii,j),1,MPIR,bufout,MPIL,ib,MPIC,ierr)
              call MPI_PACK(wsty(ii,j),1,MPIR,bufout,MPIL,ib,MPIC,ierr)
 140        continue
c
            do 150 j=js_s(n),js_e(n)
            do 150 i=is_s(n),is_e(n)
              ii=mod(i-1,NSX_M) + 1
              call MPI_PACK(ssta(ii,j),1,MPIR,bufout,MPIL,ib,MPIC,ierr)
              call MPI_PACK(sssa(ii,j),1,MPIR,bufout,MPIL,ib,MPIC,ierr)
 150        continue
 108      call MPI_BSEND(bufout,ib,MPI_PACKED,n,
     &                          MSG_MET2,MPIC,ierr)
          call M_BERR(ierr,*108)
 160    continue
c
        if(idebug.gt.2)then
          write(stderr,'(a)')
     &        ' Master:  Message MSG_MET2 sent to all slaves.'
        endif
c
c=======================================================================
c 12.  Send message type MSG_MET3 to slaves
c=======================================================================
c
        elseif(msgtype.eq.MSG_MET3)then
c  loop over slaves, clearing the input buffers on the way
          do 230 n=1,nproc
            if(lmeta(n))then
              ib = 0
              call MPI_PACK(mnextp,1,MPII,bufout,MPIL,ib,MPIC,ierr)
              do 210 j=jw_s(n),jw_e(n)
              do 210 i=iw_s(n),iw_e(n)
                ii=mod(i-1,NWX_M) + 1
                call MPI_PACK(wstx(ii,j),1,MPIR,bufout,
     &                                           MPIL,ib,MPIC,ierr)
                call MPI_PACK(wsty(ii,j),1,MPIR,bufout,
     &                                           MPIL,ib,MPIC,ierr)
 210          continue
c
              do 220 j=js_s(n),js_e(n)
              do 220 i=is_s(n),is_e(n)
                ii=mod(i-1,NSX_M) + 1
                call MPI_PACK(ssta(ii,j),1,MPIR,bufout,
     &                                           MPIL,ib,MPIC,ierr)
                call MPI_PACK(sssa(ii,j),1,MPIR,bufout,
     &                                           MPIL,ib,MPIC,ierr)
 220          continue
 109          call MPI_BSEND(bufout,ib,MPI_PACKED,n,
     &                                          MSG_MET3,MPIC,ierr)
              call M_BERR(ierr,*109)
              lmeta(n) = .false.
              nmeta = nmeta - 1
            endif
 230      continue
          lmetq = .false.
          if(nmeta.eq.0)then
            lmetp = .false.
          endif
c
       if(idebug.gt.2)then
         write(stderr,'(a)')
     &        ' Master:  Message MSG_MET3 sent to one or more slaves'
       endif
c
c=======================================================================
c 13.  Send message type MSG_MET_RDY to slaves
c=======================================================================
c
      elseif(msgtype.eq.MSG_MET_RDY)then
        do iproc=1,nproc
 110      call MPI_BSEND(bufout,1,MPI_INTEGER,
     &                          iproc,MSG_MET_RDY,MPI_COMM_WORLD,ierr)
          call M_BERR(ierr,*110)
        enddo
        if(idebug.gt.2)then
          write(stderr,'(a)')
     &         ' Master:  Message MSG_MET_RDY sent to all slaves.'
        endif
c
c=======================================================================
c  entry point to ask for tsi information (message MSG_RQ_TSI)
c  from all slaves for timestep itttsi
c=======================================================================
c
      elseif(msgtype.eq.MSG_RQ_TSI) then
        do iproc=1,nproc
 111      call MPI_BSEND(itttsi,1,MPI_INTEGER,
     &                          iproc,MSG_RQ_TSI,MPI_COMM_WORLD,ierr)
          call M_BERR(ierr,*111)
        enddo
        if(idebug.gt.2)then
          write(stderr,'(a)')
     &        ' Master:  Message MSG_RQ_TSI sent to all slaves.'
        endif
c
c=======================================================================
c  entry point to ask for part of the archive data (message MSG_RQ_ARC)
c  from selected slaves for timestep ittarch
c=======================================================================
c
      elseif(msgtype.eq.MSG_RQ_ARC) then
        ib = 0
        call MPI_PACK(ittarch,1,MPII,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(idvar  ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(jarchl ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(jarchu ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
        do i=1,nsarch
 112      call MPI_BSEND(bufout,ib,MPI_PACKED,iarcha(i),
     &                          MSG_RQ_ARC,MPI_COMM_WORLD,ierr)
          call M_BERR(ierr,*112)
        enddo
        if(idebug.gt.2)then
          write(stderr,'(2a,4i5)')
     &        ' Master:  Message MSG_RQ_ARC sent to selected slaves.',
     &        ' idvar, range = ',idvar,jarchl,jarchu,ittarch
        endif
c
c=======================================================================
c  entry point to inform the slaves that all archive data 
c  has been received and that they may clear their buffers (MSG_ARC_CLR)
c=======================================================================
c 
      elseif(msgtype.eq.MSG_ARC_CLR) then 
        do iproc=1,nproc
 113      call MPI_BSEND(ittarch,1,MPI_INTEGER,iproc,
     &                           MSG_ARC_CLR,MPI_COMM_WORLD,ierr)
          call M_BERR(ierr,*113)
        enddo
        if(idebug.gt.2)then 
          write(stderr,'(a)')
     &        ' Master:  Message MSG_ARC_CLR sent to all slaves.' 
        endif 
c
c=======================================================================
c  entry point to send snapshot request to the slaves
c=======================================================================
c
      elseif(msgtype.eq.MSG_RQ_SNAP) then
        ib = 0
        call MPI_PACK(ittsnap,1,MPII,bufout,MPIL,ib,MPIC,ierr)
        call MPI_PACK(idsnp  ,1,MPII,bufout,MPIL,ib,MPIC,ierr)
        do iproc=1,nproc
 114      call MPI_BSEND(bufout,ib,MPI_PACKED,iproc,
     &                          MSG_RQ_SNAP,MPI_COMM_WORLD,ierr)
          call M_BERR(ierr,*114)
        enddo
        if(idebug.gt.2)then
         write(stderr,'(a)')
     &        ' Master:  Message MSG_RQ_SNAP sent to all slaves.'
        endif
c
c=======================================================================
c  unrecognised request
c=======================================================================
c
      else
        write(stderr,'(a)') 
     &       ' Master:  Message type: ',msgtype,
     &       ' not yet implemented in routine m_send'
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
#include "ctmngr.h"
#include "levind.h"
#include "archive.h"
#include "snaps.h"
#include "varinfo.h"
#include "mesdta.h"
c
      integer msgtype, time
      logical lflag
      MPII = MPI_INTEGER
      MPIR = MPI_REAL
      MPIL = LENBUF*4
      MPIC = MPI_COMM_WORLD
c
c=======================================================================
c  check to see if anything to receive
c  return if no messages
c=======================================================================
c
  10  if(msgtype.ge.0)then
        call MPI_PROBE(MPI_ANY_SOURCE,msgtype,MPIC,istat,ierr)
      else
        call MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,
     &                                 MPIC,lflag,istat,ierr)
        if(.not.lflag)return
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
      msgtag = istat(MPI_TAG)
      itid = istat(MPI_SOURCE)
c
c  read name of slave and check it corresponds to the tid
c
      if(itid.lt.1.or.itid.gt.MXSLAVE)then
        write(stderr,*) ' Error: Message ',msgtag,' from tid ',itid,
     &     ' has slave number ',itid,' out of range 1 - ',MXSLAVE
c       elseif(itid.ne.itida(mes))then
c        write(6,*) ' Error: Message ',msgtag,' from slave ',mes,
c    &     ' has tid ',itid,', should be ',itida(mes)
        endif
c
c=======================================================================
c  process abort request message (MSG_BYE)
c=======================================================================
c
      if(msgtag.eq.MSG_BYE)then
        call MPI_RECV(ibufout,MPIL,MPII,itid,msgtag,MPIC,istat,ierr)
        write(stderr,*) ' Master: '
        write(stderr,*) '  Abort request received from slave: ',mes
        write(stderr,*) '  Closing down all slaves'
        call m_send(MSG_ABORT)
c
c=======================================================================
c  process position message (MSG_POSN)
c  message contains:  me, itt, itt0, iposn
c=======================================================================
c
      elseif(msgtag.eq.MSG_POSN)then
        call MPI_RECV(ibufout,MPIL,MPII,itid,msgtag,MPIC,istat,ierr)
c       call MPI_GET_COUNT(istat,MPII,lbuff,ierr)
        mes = ibufout(1)
        infoa(1,mes) = ibufout(2)
        infoa(2,mes) = ibufout(3)
        infoa(3,mes) = ibufout(4)
        if(idebug.gt.2)then
          write(stderr,901) '  Master:  Message MSG_POSN received: ',
     &    mes, infoa(1,mes),infoa(2,mes),infoa(3,mes)
 901  format(a,i5,3i10)
        endif
c
c=======================================================================
c process information message (MSG_OUTSTR)
c message contains: me, string for printing
c=======================================================================
c
      elseif(msgtag.eq.MSG_OUTSTR) then
        call MPI_RECV(bufout,MPIL,MPI_PACKED,itid,msgtag,
     &                              MPIC,istat,ierr)
        call MPI_GET_COUNT(istat,MPI_PACKED,lb,ierr)
        ib = 0
        call MPI_UNPACK(bufout,MPIL,ib,mes,1,MPII,MPIC,ierr)
        call MPI_UNPACK(bufout,MPIL,ib,outstr,80,MPI_CHARACTER
     &                                           ,MPIC,ierr)
        if(idebug.gt.2)then
          write(stderr,'(a,i5)') 
     &      '  Master:  Message MSG_OUTSTR received from slave ',mes
          write(stderr,"('  :',a,':')") outstr
        else
          write(stderr,"(' msg from ',i4,':',a)")mes,outstr
        endif
c
c=======================================================================
c  process timestep 'tsi' information (MSG_RY_TSI)
c=======================================================================
c
      elseif(msgtag.eq.MSG_RY_TSI)then
        call MPI_RECV(bufout,MPIL,MPI_PACKED,itid,msgtag,
     &                              MPIC,istat,ierr)
        call MPI_GET_COUNT(istat,MPI_PACKED,lb,ierr)
        ib = 0
        call MPI_UNPACK(bufout,MPIL,ib,mes,1,MPII,MPIC,ierr)
        call MPI_UNPACK(bufout,MPIL,ib,itttsi2,1,MPII,MPIC,ierr)
        call MPI_UNPACK(bufout,MPIL,ib,vols  ,1,MPIR,MPIC,ierr)
        call MPI_UNPACK(bufout,MPIL,ib,ektots,1,MPIR,MPIC,ierr)
        call MPI_UNPACK(bufout,MPIL,ib,dtemps,1,MPIR,MPIC,ierr)
        call MPI_UNPACK(bufout,MPIL,ib,dsalts,1,MPIR,MPIC,ierr)
c
c  checks
c
        if(.not.ltsip)
     &    write(stderr,*) ' Error: Message MSG_RY_TSI from slave ',mes,
     &                  ' sent when tsi processing inactive.' 
        if(itttsi2.ne.itttsi)then
          write(stderr,*) ' Error: Message MSG_RY_TSI from slave ',mes,
     &                    ' has wrong time flag'
          write(stderr,*) '     : itt from slave is ',itttsi2,
     &                    ' should be ', itttsi
        endif          
        if(ltsia(mes))
     &    write(stderr,*) ' Error: Message MSG_RY_TSI from slave ',mes,
     &                    ' sent twice'
        if(idebug.gt.2)then
          write(stderr,9601)'  Master:  Message MSG_RY_TSI received: ',
     &                      mes,vols,ektots,dtemps,dsalts
 9601   format(a,/,1x,i5,1p4e15.7)
        endif
c
        ltsia(mes) = .true.
        nstsi = nstsi - 1
        volume   = volume + vols
        ektot    = ektot + ektots
        dtabs(1) = dtabs(1) + dtemps
        dtabs(2) = dtabs(2) + dsalts
        if(idebug.gt.2)then
          write(stderr,9602) 
     &         '  Master:  Message MSG_RY_TSI. volume, ektot, dtabs = ',
     &          volume, ektot, dtabs(1), dtabs(2)
 9602   format(a,/,1x,1p4e15.7)
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
          total_time=rtc()*6.6e-9-start_time
          write(stderr,9603) itttsi,timetsi,ektot,dtabs(1),dtabs(2),
     &                       total_time
9603      format(1x,'ts=',i7,1x,a32,' ke=', 1pe13.6,
     &    ' dtemp=',1pe13.6,' dsalt=',1pe13.6,' time=',0pf10.4)
#else
          itime = MPI_WTIME()-start_time
          write(stderr,9603) itttsi,timetsi,ektot,dtabs(1),dtabs(2),
     &                       itime
9603      format(1x,'ts=',i7,1x,a32,' ke=', 1pe13.6,
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
        call MPI_RECV(bufout,MPIL,MPIR,itid,msgtag,MPIC,istat,ierr)
        call MPI_GET_COUNT(istat,MPIR,lbuff,ierr)
        lbuff = nint(bufout(1))
        if(lbuff.gt.LENBUF)then
           write(stderr,*) ' error in m_recv when receiving ',
     &                      ' MSG_RY_ARC'
           write(stderr,*) ' buffer overflow'
        endif
        ittsav  = nint(bufout(2))
        idvar2  = nint(bufout(3))
        lbuff   = 3
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
 55    imast = nint(bufout(lbuff+1))
c
c check for "end of message" marker
c
       if(imast.gt.0) then
         lbuff = lbuff+2
         jmast = nint(bufout(lbuff))
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
            lbuff = lbuff+1
            rest2d(imast,jmast) = bufout(lbuff)
          elseif(vartgrd(idvar2)) then
c
c convert jmast to sub-volume index
c
            jzon = jmast - ((jmast-1)/JSUB_M)*JSUB_M
            do 56 k=1,kmt(imast,jmast)
              lbuff = lbuff+1
              rest3d(k,imast,jzon) = bufout(lbuff)
  56        continue
          else
c
c convert jmast to sub-volume index
c
           jzon = jmast - ((jmast-1)/JSUB_M)*JSUB_M 
           do 57 k=1,kmu(imast,jmast)
             lbuff = lbuff+1
             rest3d(k,imast,jzon) = bufout(lbuff)
  57       continue
         endif   
         goto 55
       endif
c
       if(nseapt.ge.nseachk) call archwr       
c
       if(idebug.gt.2)then
        write(stderr,'(a,i5,a,4i5)') 
     &       '  Master:  Message MSG_RY_ARC from slave ',mes,
     &       '  idvar, itt, nseapt, nseachk = ',
     &          idvar,ittsav,nseapt,nseachk
       endif
c
c=======================================================================
c  process snapshot replies (MSG_RY_SNAP)
c=======================================================================
c
      elseif(msgtag.eq.MSG_RY_SNAP)then
        call MPI_RECV(bufout,MPIL,MPIR,itid,msgtag,MPIC
     &                                              ,istat,ierr)
        call MPI_GET_COUNT(istat,MPIR,lbuff,ierr)
        lbuff = nint(bufout(1))
        if(lbuff.gt.LENBUF)then
           write(stderr,*) ' error in m_recv when receiving ',
     &                      ' MSG_RY_SNAP'
           write(stderr,*) ' buffer overflow'
        endif
        ittsnap2  = nint(bufout(2))
        idsnp2    = nint(bufout(3))
        lbuff     = 3
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
 65      imast = nint(bufout(lbuff+1))
         if(imast.gt.0)then
           jmast = nint(bufout(lbuff+2))
c
c check for corrupted message
c
          if(imast.gt.IMT_M.or.jmast.gt.JMT_M) then
            write(stderr,*) ' Invalid master indices: ',imast,jmast
            write(stderr,*) ' currently unpacking ',idsnp,'th snapshot'
            call m_send(MSG_ABORT)
          endif
          snap2d(imast,jmast) = bufout(lbuff+3)
          lbuff = lbuff+3
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
        write(stderr,*) '  Master:  Message MSG_RY_SNAP from slave ',
     &                    mes,'  idsnp, itt, nssnap = ',
     &                    idsnp2,ittsnap2,nssnap         
        endif
c
c=======================================================================
c  receive met requests (MSG_CLR_MET)
c=======================================================================
c
      elseif(msgtag.eq.MSG_CLR_MET)then
        call MPI_RECV(ibufout,MPIL,MPII,itid,msgtag,MPIC,istat,ierr)
        mes   = ibufout(1)
        mplus = ibufout(2)
        if(lmeta(mes))then
          write(stderr,*)' Master: Message MSG_CLR_MET received from ',
     &           'slave ',mes,' before previous request answered'
        endif
        if(mplus.ne.mnextp)then
          write(stderr,*)' Master: Message MSG_CLR_MET requests wrong ',
     &           'month.  Request for:',mplus,',   should be:',mnextp
        endif
c
c  set flags so that main loop sends data to requesting slave
c
        lmetq      = .true.
        lmeta(mes) = .true.
c
        if(idebug.gt.2)then
          write(stderr,'(a,i5)') 
     &       '  Master:  Message MSG_CLR_MET received from slave ', mes 
        endif
c
c=======================================================================
c  error          
c=======================================================================
c
      else
        write(stderr,*) ' Master (m_recv):  Message type: ',msgtag,
     &                  ' received from slave',mes
        write(stderr,*) ' This message type is not yet implemented'
      endif
c
c=======================================================================
c  flush stderr to ensure debug messages are printed
c=======================================================================
c
      if(idebug.gt.2)then
        call flush(stderr)
      endif
c
c=======================================================================
c  PART 2.
c  if OK to send, send off any flux messages now ready
c  to go
c=======================================================================
c
      if(.not.lsend)return
c
      if(lmetp .and. lmetq)then
        call m_send(MSG_MET3)
      endif
c
c=======================================================================
c  flush stderr to ensure debug messages are received
c=======================================================================
c
      if(idebug.gt.2)then
          call flush(stderr)
      endif
c
      if(msgtype.ne.msgtag)goto 10
      return
      end


      subroutine M_BERR(ierr,*)
c
c  MPI master subroutine to check and process return code from MPI_BSEND
c
#include "def_slave.h"
#include "param.h"
#include "mesdta.h"
c
      character*( MPI_MAX_ERROR_STRING ) string
      double precision time0,time1
      save in,time0
      data in/0/
c
      if(ierr.eq.0)then
        lsend = .true.
        in = 0
        return
      endif
c
c  print error message on first entry
c
      if(in.eq.0)then
        lsend = .false.
        in    = 1
        time0 = MPI_WTIME()
        call MPI_ERROR_CLASS(ierr,ierrc)
        call MPI_ERROR_STRING(ierr,string,len,ierr2)
        write(stderr,10)ierr,ierrc,time0,string(1:len)
        call flush(stderr)
        call m_recv(MSG_ANY)
        return 1
 10   format(' Master: Error with MPI_BSEND ',/,
     &      '   ierr, ierrc, time = ',2i10,f16.4,/,
     &      '   error string      = ',a)
      endif
c
c  return again ?
c
      in    = in+1
      time1 = MPI_WTIME()
      if(time1-time0.gt.10d0)then
        call m_recv(MSG_ANY)
        return 1
      endif
c
c  abort
c
      call MPI_ERROR_CLASS(ierr,ierrc)
      call MPI_ERROR_STRING(ierr,string,len,ierr2)
      write(stderr,11)ierr,ierrc,time1,string(1:len)
 11   format(' Master: Continuing error with MPI_BSEND',/,
     &      '   ierr, ierrc, time = ',2i10,f16.4,/,
     &      '   error string      = ',a,
     &      ' Program aborting ...')
c
c  exit
c
      call flush(stderr)
      call MPI_ABORT(MPI_COMM_WORLD,-999)
      stop
      end
        
