      subroutine initmess
c
c @(#) SCCS module: s_msg.F  version: 1.3
c     Creation date: 03/17/95
c
c=======================================================================
c  slave subroutine to initialise message passsing logic
c=======================================================================
c
#include "def_slave.h"
#include "param.h"
#include "mesdta.h"
c
c  set logical values
c
      lbrkpt = .false.
      ltsip  = .false.
      ltsiw  = .false.
      lsnapp = .false.
      larchp = .false.
c
c  set flags to indicate that all boundary data is in place at the
c  start of the first timestep
c
      nfalse(1)=0
      nfalse(2)=0
      iclnclr = 0
      lbcln   = .false.
      lbtrp   = .false.
      itrpclr = 0
c
      return
      end


      subroutine schkpnt(iposn0)
c
c=======================================================================
c  slave subroutine to report to master and wait until 
c  message to proceed has arrived
c=======================================================================
c
c  input:  iposn0 - position 
c
#include "param.h"
#include "mesdta.h"
c
c  send position message MSG_POSN
c
      iposn = iposn0
      call s_send(MSG_POSN)
c
c  process messages until continue (message MSG_CONT) received 
c  and flag set
c
      lbrkpt = .false.
 10   call s_recv(MSG_ANY)
      if(.not.lbrkpt)goto 10 
      return
      end

      subroutine sendposn(iposn0)
c
c=======================================================================
c  slave subroutine to report position to master 
c=======================================================================
c
c  input:  iposn0 - position 
c
#include "param.h"
#include "mesdta.h"
c
c  send position message MSG_POSN
c
      iposn = iposn0
      call s_send(MSG_POSN)
c
c  exit
c
      return
      end

      subroutine initslvc(n)
c
c=======================================================================
c  initialise slave flags ready to receive 'ready to send' messages
c=======================================================================
c
#include "param.h"
#include "mesdta.h"
c
      if(n.eq.1) then
        lbcln = .false.
        iclnclr = 0
      else
        lbtrp = .false.
        itrpclr = 0
      endif
      return
      end

      subroutine initslvd(n)
c
c=======================================================================
c  initialise slave message passing logic ready to receive new data
c=======================================================================
c
c     iota = true for all but outer halo points at all times
c          = false for outer halo points before data arrives
c
#include "param.h"
#include "ctmngr.h"
#include "mesdta.h"
#include "timelv.h"
c
      do 10 j = 1, JMT_S
      do 10 i = 1, IMT_S
        iota(i,j,n) = .true.
   10 continue
c
      do 20 k = 1, nouter
        i = iout(k)
        j = jout(k)
        iota(i,j,n) = .false.
   20 continue
      nfalse(n) = nouter
c
      if(n.eq.1) then
        do 30 na = 1, nactive
         bclnflg(na)=.false.
   30   continue
        itt3 = itt
        np3  = np
        if(mxpas2)then
          ittp3 = 2
        else
          ittp3 = 1
        endif
      else
        do 40 na = 1, nactive 
          btrpflg(na)=.false. 
   40   continue 
        itt2   = itt
        itbt2  = itbt
        itbtp2 = itbtp
        np02   = np0
      endif
c
      return
      end


      subroutine msgbye
c
c=======================================================================
c subroutine to send the MSG_BYE message to the master process
c This subroutine is a copy of the entry in the block-if in
c s_send and should only be used to send the message from within
c s_send in order to avoid recursion. All other instances should
c use the "call s_send(MSG_BYE)" statement
c=======================================================================
c
#include "param.h"
#include "mesdta.h"
#include "switch.h"
c
c Note: pvmfinitsend will clear the buffer of any partially packed 
c messages.
c
      call pvmfinitsend( PVMDEFAULT, info )
      call pvmfpack( INTEGER_PVM,    me, 1, 1, info )
      call pvmfsend( mtid, MSG_BYE, info )
      if(idebug.gt.3)then
        write(ioslave,*)' Slave - message MSG_BYE sent. '
      endif
      call schkpnt(1)
      call closedown(info)
      stop
      end


      subroutine s_send(msgtype)
c
c=======================================================================
c  generic message handling subroutine for moma.pvm.  This
c  subroutine is used by the slave programs to pack 
c  and send messages.
c=======================================================================
c
#include "def_slave.h"
#include "param.h"
#include "scalar.h"
#include "slabs.h"
#include "frees.h"
#include "switch.h"
#include "iounit.h"
#include "cdiag.h"
#include "chmix.h"
#include "ctmngr.h"
#include "timelv.h"
#include "levind.h"
#include "archive.h"
#include "snaps.h"
#include "varinfo.h"
#include "mesdta.h"
      integer msgtype
c
c=======================================================================
c  Type MSG_BYE: inform the master that the slave is aborting
c=======================================================================
c
      if(msgtype.eq.MSG_BYE)then
       call pvmfinitsend( PVMDEFAULT, info )
       call pvmfpack( INTEGER_PVM,    me, 1, 1, info )
       call pvmfsend( mtid, msgtype, info )
       if(idebug.gt.3)then
         write(ioslave,*)' Slave - message MSG_BYE sent. '
       endif
       call schkpnt(1)
       call closedown(info)
       stop
c
c=======================================================================
c  Type MSG_POSN: send status information message
c=======================================================================
c
      elseif(msgtype.eq.MSG_POSN)then
       call pvmfinitsend( PVMDEFAULT, info )
       call pvmfpack( INTEGER_PVM,    me, 1, 1, info )
       call pvmfpack( INTEGER_PVM,    iposn, 1, 1, info )
       call pvmfpack( INTEGER_PVM,    itt, 1, 1, info )
       call pvmfpack( INTEGER_PVM,    itbt, 1, 1, info )
       call pvmfsend( mtid, msgtype, info )
       if(idebug.gt.3)then
         write(ioslave,*)' Slave - message MSG_POSN sent.',
     &                   ' iposn=',iposn
       endif
c
c=======================================================================
c Type MSG_OUTSTR: send string 'outstr' for printing by master
c=======================================================================
c
      elseif(msgtype.eq.MSG_OUTSTR) then
       call pvmfinitsend( PVMDEFAULT, info )
       call pvmfpack( INTEGER_PVM,    me, 1, 1, info )
       call pvmfpack( STRING,   outstr, 80, 1, info )
       call pvmfsend( mtid, msgtype, info )
       if(idebug.gt.3)then
         write(ioslave,*)' Slave - message MSG_OUTSTR sent.',
     &                   ' outstr =',outstr
       endif
c
c=======================================================================
c Type MSG_RY_TSI: pack and send time step information to master
c=======================================================================
c 
      elseif(msgtype.eq.MSG_RY_TSI) then
       call pvmfinitsend( PVMDEFAULT, info )
       call pvmfpack( INTEGER_PVM,     me, 1, 1, info )
       call pvmfpack( INTEGER_PVM,    itt, 1, 1, info )
       call pvmfpack( REAL_PVM,    volume, 1, 1, info )
       call pvmfpack( REAL_PVM,     ektot, 1, 1, info )
       call pvmfpack( REAL_PVM,  dtabs(1), 2, 1, info )
       call pvmfpack( REAL_PVM,  dtabs(1), 2, 1, info )
       call pvmfsend( mtid, msgtype, info )
c unset tsi flag (in case next master request arrives early)
       prntsn = .false.
       if(idebug.gt.3)then
         write(ioslave,*)' Slave - message MSG_RY_TSI sent'
       endif
c
c=======================================================================
c Type MSG_RY_ARC: pack and send partial archive reply to master
c=======================================================================
c 
      elseif(msgtype.eq.MSG_RY_ARC) then
c
c convert jarchu and jarchl to "slave indices" and limit to slave bounds
c
        jarchl = max(1,jarchl-jswm1)
        jarchu = min(jt,jarchu-jswm1)
        nchk = 0
c
c pack and send requested data
c
         do 50 j=jarchl,jarchu
         do 50 i=2,itm1
          if(kpn(i,j).eq.me) then
c
c  initialise buffer
c
            if(nchk.eq.0 .and.(vartgrd(idvar).or.kmu(i,j).ne.0))then
              call pvmfinitsend( PVMDEFAULT, info )
              call pvmfpack(INTEGER_PVM,    me, 1, 1, info )
              call pvmfpack(INTEGER_PVM,ittsav, 1, 1, info)
              call pvmfpack(INTEGER_PVM, idvar, 1, 1, info)
            endif
c
c convert to master indices and pack values
c
            imast = mod(iswm1 + i-1,IMT_M) + 1
            jmast = jswm1 + j
c
c Now pack requested variable(s) for this point
c 
c
            if(idvar.le.NUM2D.and.vartgrd(idvar)) then
              nchk =nchk + 1
              call pvmfpack(INTEGER_PVM, imast, 1, 1, info)
              call pvmfpack(INTEGER_PVM, jmast, 1, 1, info)
              call pvmfpack(REAL_PVM, hsave(i,j,idvar), 1, 1, info)
            elseif(idvar.le.NUM2D.and.kmu(i,j).ne.0) then
              nchk =nchk + 1
              call pvmfpack(INTEGER_PVM, imast, 1, 1, info)
              call pvmfpack(INTEGER_PVM, jmast, 1, 1, info)
              call pvmfpack(REAL_PVM, hsave(i,j,idvar), 1, 1, info)
            elseif(vartgrd(idvar)) then
              nchk =nchk + 1
              call pvmfpack(INTEGER_PVM, imast, 1, 1, info)
              call pvmfpack(INTEGER_PVM, jmast, 1, 1, info)
              call pvmfpack(REAL_PVM, tsave(1,i,j,idvar-NUM2D), 
     &                        kmt(i,j), 1, info)
            elseif(kmu(i,j).ne.0) then
              nchk =nchk + 1
              call pvmfpack(INTEGER_PVM, imast, 1, 1, info)
              call pvmfpack(INTEGER_PVM, jmast, 1, 1, info)
              call pvmfpack(REAL_PVM, tsave(1,i,j,idvar-NUM2D), 
     &                        kmu(i,j), 1, info)
            endif
          endif
  50    continue
c
c  send buffer only if it has data in it
c
        if(nchk.ne.0)then
c
c Pack a negative integer to indicate "end of message"
c
          call pvmfpack(INTEGER_PVM,  -999, 1, 1, info)
          call pvmfsend( mtid, MSG_RY_ARC, info )
          if(idebug.gt.3) then
            write(ioslave,*)' Slave - message MSG_RY_ARC sent',
     &                 ' for idvar:',idvar,'  nsea chk = ',nchk
          endif
        else
          if(idebug.gt.3) then
            write(ioslave,*)' Slave - message MSG_RY_ARC not sent.',
     &          ' for idvar:',idvar, '  Slave had no data to send.'
           endif
        endif
c
c  clear the "partial out-standing flag"
c
        larchq = .false.
c
c=======================================================================
c Type MSG_RY_SNAP: pack and send partial snapshot reply to master
c=======================================================================
c
      elseif(msgtype.eq.MSG_RY_SNAP) then
c
c
c pack and send requested data
c
         idsvar=ispvar(idsnp)
         kslev =isplev(idsnp)
         nchk = 0
c
         do 60 j=1,jt
         do 60 i=2,itm1
          if(kpn(i,j).eq.me.and.(kslev.le.kmu(i,j) .or. 
     &         (vartgrd(idsvar).and.kslev.le.kmt(i,j)) )  ) then
c
c  initialise buffer
c
            if(nchk.eq.0)then
              nchk = 1
              call pvmfinitsend( PVMDEFAULT, info )
              call pvmfpack(INTEGER_PVM,    me, 1, 1, info )
              call pvmfpack(INTEGER_PVM,ittsnap, 1, 1, info)
              call pvmfpack(INTEGER_PVM, idsnp, 1, 1, info)
            endif
c
c convert to master indices and pack values
c
            imast = mod(iswm1 + i-1,IMT_M) + 1
            jmast = jswm1 + j
c
c Now pack requested variable(s) for this point
c
            call pvmfpack(INTEGER_PVM, imast, 1, 1, info)
            call pvmfpack(INTEGER_PVM, jmast, 1, 1, info)
            call pvmfpack(REAL_PVM, snap2d(i,j,idsnp), 1, 1, info)
          endif
  60     continue
         if(nchk.ne.0)then
c
c Pack a negative integer to indicate "end of message"
c
           call pvmfpack(INTEGER_PVM,  -999, 1, 1, info)
           call pvmfsend( mtid, MSG_RY_SNAP, info )
           if(idebug.gt.3) then
             write(ioslave,*)' Slave . snapshot reply sent to',
     &                       ' master.  idsnp = ', idsnp
           endif
         else
           if(idebug.gt.3) then
            write(ioslave,*)' Slave - message MSG_RY_SNAP not sent.',
     &          ' for idvar:',idvar, '  Slave had no data to send.'
           endif
         endif
c
c  clear the "partial out-standing flag" plus "snapshot in progress" 
c  and snapshot needed this timestep" flags if this is the last field
c
         lsnapq = .false.
         if(idsnp.eq.NSNAPS)then
c unset snapts flag (in case next master request arrives early)
           snapts = .false.
           lsnapp = .false.
         endif
c
c
c=======================================================================
c Type MSG_BOUND1: pack and send 3-d boundary information to other
c slaves
c=======================================================================
c
      elseif (msgtype .eq. MSG_BOUND1) then
       moffset=0
       do 20 knb = nactive, 1, -1
         mes=nactme(knb)
         call pvmfinitsend( PVMDEFAULT, info)
#ifdef pvm_buffer
         bufout(2) = me
         bufout(3) = itt3
         bufout(4) = ittp3
         bufout(5) = nsend(mes)
         lbuff = 5
#else
         call pvmfpack( INTEGER_PVM, me, 1, 1, info)
         call pvmfpack( INTEGER_PVM, itt3, 1, 1, info)
         call pvmfpack( INTEGER_PVM, ittp3, 1, 1, info)
         call pvmfpack( INTEGER_PVM, nsend(mes), 1, 1, info)
#endif
         do 10 m = moffset+1, moffset+nsend(mes)
          if (isd(m,3) .ne. mes) then
            write(ioslave,*) ' s_send.  error while packing',
     &                       ' MSG_BOUND1:'
            write(ioslave,*) ' knb,mes,m,moffset,nsend(mes) = ',
     &                         knb,mes,m,moffset,nsend(mes)
            write(ioslave,*) ' isd(m,1),isd(m,2),isd(m,3) = ',
     &                         isd(m,1),isd(m,2),isd(m,3)
           call msgbye
          endif
          i = mod(isd(m,1)-1 + iswm1, IMT_M) + 1
          j = isd(m,2) + jswm1
#ifdef pvm_buffer
          lbuff = lbuff + 2
          if(lbuff+kmt(isd(m,1), isd(m,2))*2
     &            +kmu(isd(m,1), isd(m,2))*2 .gt.LENBUF)then
            write(ioslave,*) ' error in s_send while packing',
     &                       ' MSG_BOUND1:'
            write(ioslave,*) ' buffer overflow'
            write(ioslave,*) ' m,nsend(mes)= ', m,nsend(mes)
            call msgbye
          endif
          bufout(lbuff-1) = i
          bufout(lbuff  ) = j
          do 51 k=1, kmt(isd(m,1), isd(m,2))
            lbuff = lbuff+1
            bufout(lbuff) = t(k,isd(m,1),isd(m,2),1,np)
  51      continue
          do 52 k=1, kmt(isd(m,1), isd(m,2))
            lbuff = lbuff+1
            bufout(lbuff) = t(k,isd(m,1),isd(m,2),2,np)
  52      continue
          do 53 k=1, kmu(isd(m,1), isd(m,2))
            lbuff = lbuff+1
            bufout(lbuff) = u(k,isd(m,1),isd(m,2),np)
  53      continue
          do 54 k=1, kmu(isd(m,1), isd(m,2))
            lbuff = lbuff+1
            bufout(lbuff) = v(k,isd(m,1),isd(m,2),np)
  54      continue
#else
          call pvmfpack(INTEGER_PVM, i, 1, 1, info)
          call pvmfpack(INTEGER_PVM, j, 1, 1, info)
          call pvmfpack(REAL_PVM, t(1,isd(m,1),isd(m,2),1,np),
     &                  kmt(isd(m,1), isd(m,2)), 1, info)
          call pvmfpack(REAL_PVM, t(1,isd(m,1),isd(m,2),2,np),
     &                  kmt(isd(m,1),isd(m,2)), 1, info)
          call pvmfpack (REAL_PVM, u(1,isd(m,1),isd(m,2),np),
     &                  kmu(isd(m,1), isd(m,2)), 1, info)
          call pvmfpack (REAL_PVM, v(1,isd(m,1),isd(m,2),np),
     &                  kmu(isd(m,1), isd(m,2)), 1, info)
#endif
   10    continue
#ifdef pvm_buffer
         bufout(1) = lbuff
         call pvmfpack (REAL_PVM, bufout, lbuff, 1, info)
#endif
         call pvmfsend( itida(mes), msgtype, info)
       moffset=moffset + nsend(mes)
   20  continue
       if(idebug.gt.4)then
         write(ioslave,*)' Slave - message MSG_BOUND1 sent.',
     &                   ' iposn=',iposn
       endif
c
c=======================================================================
c Type MSG_BOUND2: pack and send 2-d boundary information to other
c     slaves.  Free surface timestep is incremented before sending
c=======================================================================
c
      elseif (msgtype .eq. MSG_BOUND2) then
       moffset = 0
       do 40 knb = nactive, 1, -1
        mes = nactme(knb)
        call pvmfinitsend( PVMDEFAULT, info)
#ifdef pvm_buffer
        bufout(2) = me
        bufout(3) = itt2
        bufout(4) = itbt2
        bufout(5) = itbtp2
        bufout(6) = nsend(mes)
        lbuff = 6
#else
        call pvmfpack( INTEGER_PVM, me, 1, 1, info)
        call pvmfpack( INTEGER_PVM, itt2, 1, 1, info)
        call pvmfpack( INTEGER_PVM, itbt2, 1, 1, info)
        call pvmfpack( INTEGER_PVM, itbtp2, 1, 1, info)
        call pvmfpack( INTEGER_PVM, nsend(mes), 1, 1, info)
#endif
        do 30 m = moffset + 1, moffset+nsend(mes)
          if (isd(m,3) .ne. mes) then
           write(ioslave,*) ' s_send: error while packing MSG_BOUND2'
           write(ioslave,*) ' isd,m,mes = ', 
     &                      isd(m,1),isd(m,2),isd(m,3),m,mes
           call msgbye
          endif
          i = mod(isd(m,1) - 1 + iswm1, IMT_M) + 1
          j = isd(m,2) + jswm1
#ifdef pvm_buffer
          lbuff = lbuff + 5
          if(lbuff.gt.LENBUF)then
            write(ioslave,*) ' s_send: error while packing MSG_BOUND2'
            write(ioslave,*) ' buffer overflow'
            write(ioslave,*) ' m,nsend(mes)= ', m,nsend(mes)
            call msgbye
          endif
          bufout(lbuff-4) = i
          bufout(lbuff-3) = j
          bufout(lbuff-2) = h0(isd(m,1),isd(m,2),np0)
          bufout(lbuff-1) = u0(isd(m,1),isd(m,2),np0)
          bufout(lbuff  ) = v0(isd(m,1),isd(m,2),np0)
#else
          call pvmfpack (INTEGER_PVM, i, 1, 1, info)
          call pvmfpack (INTEGER_PVM, j, 1, 1, info)
          call pvmfpack (REAL_PVM, h0(isd(m,1),isd(m,2),np0), 1, 1, info)
          call pvmfpack (REAL_PVM, u0(isd(m,1),isd(m,2),np0), 1, 1, info)
          call pvmfpack (REAL_PVM, v0(isd(m,1),isd(m,2),np0), 1, 1, info)
#endif
   30   continue
#ifdef pvm_buffer
        bufout(1) = lbuff
        call pvmfpack (REAL_PVM, bufout, lbuff, 1, info)
#endif
        call pvmfsend( itida(mes), msgtype, info)
       moffset = moffset + nsend(mes)
   40 continue
      if(idebug.gt.4)then
        write(ioslave,*)' Slave - message MSG_BOUND2 sent. iposn=',
     &                    iposn
      endif
c
c=======================================================================
c Type MSG_CLR_BD1: inform active neighbours that you are ready for
c      for the next batch of baroclinic boundary information
c=======================================================================
c
      elseif (msgtype .eq. MSG_CLR_BD1) then
        call pvmfinitsend( PVMDEFAULT, info)
        call pvmfpack( INTEGER_PVM, me, 1, 1, info)
        call pvmfmcast( nactive, nactid, MSG_CLR_BD1, info )
        if(idebug.gt.4)then
          write(ioslave,*)
     &    ' Slave. Message MSG_CLR_BD1 sent to active neighbours'
        endif
c 
c=======================================================================
c Type MSG_CLR_BD2: inform active neighbours that you are ready for 
c      for the next batch of barotropic boundary information
c=======================================================================
c
      elseif (msgtype .eq. MSG_CLR_BD2) then 
        call pvmfinitsend( PVMDEFAULT, info) 
        call pvmfpack( INTEGER_PVM, me, 1, 1, info) 
        call pvmfmcast( nactive, nactid, MSG_CLR_BD2, info ) 
        if(idebug.gt.4)then 
         write(ioslave,*)
     &    ' Slave. Message MSG_CLR_BD2 sent to active neighbours'
        endif
c
c=======================================================================
c Unrecognised message type:
c=======================================================================
c
      else
c
       write(outstr,901) msgtype
 901  format(' Slave - trying to send message type: ',i5,
     &                 '.    Not yet implemented.')
       write(ioslave,*) outstr
       call pvmfinitsend( PVMDEFAULT, info )
       call pvmfpack( INTEGER_PVM,     me, 1, 1, info )
       call pvmfpack( STRING,   outstr, 80, 1, info )
       call pvmfsend( mtid, MSG_OUTSTR, info )
       if(idebug.gt.3)then
         write(ioslave,*)' Slave - message MSG_OUTSTR indicating',
     &                   ' error sent'
       endif
      endif
      return
      end


      subroutine s_recv(msgtype)
c
c=======================================================================
c  generic message handling subroutine for moma.pvm.  This
c  subroutine is used by the slave programs to wait for, 
c  receive and unpack messages from the master.
c=======================================================================
c
c  if msgtype .ge. zero:
c      1.  wait for required message and exit aftr processing
c  otherwise:    
c      1.  check for, receive and unpack any waiting messages,
c          and either set flags or respond.
c      2.  check the message passing flags and send off messages
c          which can now be sent
c
#include "def_slave.h"
#include "param.h"
#include "scalar.h"
#include "frees.h"
#include "slabs.h"
#include "switch.h"
#include "iounit.h"
#include "ctmngr.h"
#include "timelv.h"
#include "chmix.h"
#include "cvmix.h"
#include "grdvar.h"
#include "coord.h"
#include "levind.h"
#include "archive.h"
#include "snaps.h"
#include "varinfo.h"
#include "mesdta.h"
c
      integer msgtype
c
c=======================================================================
c   PART 1.  blocking or non-blocking receive.
c=======================================================================
c 
c
c=======================================================================
c  check for abort requests from the Master
c=======================================================================
c
        call pvmfprobe(mtid, MSG_ABORT, ibabort)
        if (ibabort.gt.0) then
           if(idebug.gt.3)then
             write(ioslave,*)' Slave - message MSG_ABORT received'
           endif
           call sendposn(1)
           call closedown(info)
           stop
        endif
c
  10  if(msgtype.ge.0)then
        call pvmfrecv  (  -1, msgtype, ibuff)
      else
        call pvmfnrecv (  mtid, MSG_ANY, ibuff)
      endif
c
c=======================================================================
c  processes messages 
c=======================================================================
c
      if(ibuff.ne.0)then
        call pvmfbufinfo(ibuff, ibytes, msgtag, itid, info)
#ifdef cray-t3d
        call pvmfgetpe(itid,itidpe)
        itid = itidpe
#endif
c
c=======================================================================
c  0.  Receive abort and close-down message
c=======================================================================
c
        if(msgtag.eq.MSG_ABORT) then
          if(idebug.gt.3)then
            write(ioslave,*)' Slave - message MSG_ABORT received'
          endif
          call sendposn(1)
          call closedown(info)
          stop
c
c=======================================================================
c  1.  Receive message type MSG_CONT
c      wait for permission to continue and set flag
c
c=======================================================================
c
        elseif(msgtag.eq.MSG_CONT) then
          lbrkpt = .true.
          if(idebug.gt.3)then
            write(ioslave,*)' Slave - message MSG_CONT received'
          endif
c  
c=======================================================================
c  2.  Receive message type MSG_IDS
c      receive and unpack tid information from the master
c=======================================================================
c
        elseif(msgtag.eq.MSG_IDS) then
          call pvmfunpack( INTEGER_PVM, nproc, 1, 1, info )
          call pvmfunpack( INTEGER_PVM, itida, nproc, 1, info)
          if(idebug.gt.3)then
            write(ioslave,*)' Slave - message MSG_IDS received'
          endif
c
c=======================================================================
c  3.  Receive message type MSG_CONTROL
c      receive and unpack control information from the master
c=======================================================================
c
        elseif(msgtag.eq.MSG_CONTROL) then
          call pvmfunpack( INTEGER_PVM,    init, 1, 1, info )
          call pvmfunpack( INTEGER_PVM,  restrt, 1, 1, info )
          call pvmfunpack( INTEGER_PVM,      eb, 1, 1, info )
          call pvmfunpack( REAL_PVM,     days, 1, 1, info )
          call pvmfunpack( REAL_PVM,      tsi, 1, 1, info )
          call pvmfunpack( REAL_PVM,     acor, 1, 1, info )
          call pvmfunpack( REAL_PVM,   dgnstc, 1, 1, info )
          call pvmfunpack( REAL_PVM,    snapd, 1, 1, info )
          call pvmfunpack( REAL_PVM,    archd, 1, 1, info )
          call pvmfunpack( REAL_PVM,     dtts, 1, 1, info )
          call pvmfunpack( REAL_PVM,     dtuv, 1, 1, info )
          call pvmfunpack( REAL_PVM,     dtbt, 1, 1, info )
          call pvmfunpack( REAL_PVM,       am, 1, 1, info )
          call pvmfunpack( REAL_PVM,       ah, 1, 1, info )
          call pvmfunpack( REAL_PVM,     fkpm, 1, 1, info )
          call pvmfunpack( REAL_PVM,     fkph, 1, 1, info )
          call pvmfunpack( REAL_PVM,    cdbot, 1, 1, info )
          call pvmfunpack( REAL_PVM,    dxdeg,    1, 1, info )
          call pvmfunpack( REAL_PVM,    dydeg,    1, 1, info )
          call pvmfunpack( REAL_PVM,       dz,   KM, 1, info )
#ifdef de_checkbd
          call pvmfunpack( REAL_PVM,   dchkbd, 1, 1, info )
#endif
          call pvmfunpack( INTEGER_PVM,  nmix, 1, 1, info )
          call pvmfunpack( INTEGER_PVM,  ncon, 1, 1, info )
          call pvmfunpack( INTEGER_PVM,  ntbt, 1, 1, info )
          call pvmfunpack( INTEGER_PVM,idebug, 1, 1, info )
          call pvmfunpack( INTEGER_PVM,iorest, 1, 1, info )
          call pvmfunpack( INTEGER_PVM,ispvar, NSNAPS, 1, info )
          call pvmfunpack( INTEGER_PVM,isplev, NSNAPS, 1, info )
          call pvmfunpack( STRING, fnrest,80, 1, info )
          if(idebug.gt.3)then
            write(ioslave,*)' Slave - message MSG_CONTROL received'
          endif
c
c=======================================================================
c  4.  Receive message type MSG_TOPOG
c=======================================================================
c
      elseif(msgtag.eq.MSG_TOPOG) then
c-----------------------------------------------------------------------
c receive and unpack kmt and grid information from the master
c-----------------------------------------------------------------------
          call pvmfunpack( INTEGER_PVM,    it,   1, 1, info )
          call pvmfunpack( INTEGER_PVM,    jt,   1, 1, info )
          call pvmfunpack( INTEGER_PVM, icycl,   1, 1, info )
          call pvmfunpack( INTEGER_PVM, iswm1,   1, 1, info )
          call pvmfunpack( INTEGER_PVM, jswm1,   1, 1, info )
          call pvmfunpack( REAL_PVM,    stlon,   1, 1, info )
          call pvmfunpack( REAL_PVM,    stlat,   1, 1, info )
          do 20 j=1,jt
          do 20 i=1,it
           call pvmfunpack( INTEGER_PVM, kmt(i,j),  1, 1, info )
 20       continue
          do 30 j=1,jt
          do 30 i=1,it
           call pvmfunpack( INTEGER_PVM, kmu(i,j),  1, 1, info )
 30       continue
          do 40 j=1,jt
          do 40 i=1,it
           call pvmfunpack( INTEGER_PVM, kpn(i,j),  1, 1, info )
 40       continue
          if(idebug.gt.3)then
            write(ioslave,*)' Slave - message MSG_TOPOG received'
          endif
c
c=======================================================================
c  6.  Receive message type MSG_RQ_TSI
c      timestep 'tsi' information request from master
c=======================================================================
c
        elseif(msgtag.eq.MSG_RQ_TSI)then
          call pvmfunpack( INTEGER_PVM,itttsi,    1, 1, info )
          ltsip = .true.
          if(idebug.gt.3)then
            write(ioslave,*)' Slave - message MSG_RQ_TSI received'
          endif
c
c=======================================================================
c  7.  Receive message type MSG_RQ_ARC
c      timestep 'tsi' information request from master
c=======================================================================
c
        elseif(msgtag.eq.MSG_RQ_ARC)then
          larchq = .true.
          call pvmfunpack( INTEGER_PVM, ittarch, 1, 1, info )
          call pvmfunpack( INTEGER_PVM, idvar  , 1, 1, info )
          call pvmfunpack( INTEGER_PVM, jarchl , 1, 1, info )
          call pvmfunpack( INTEGER_PVM, jarchu , 1, 1, info )
          if(idebug.gt.3)then
            write(ioslave,*)' Slave - message MSG_RQ_ARC received,',
     &                      ' idvar = ',idvar
            write(outstr,900) idvar
            call s_send(MSG_OUTSTR)
  900  format(' Slave - message MSG_RQ_ARC received, idvar = ',I10)
          endif
c 
c=======================================================================
c  8.  Receive message type MSG_ARC_CLR 
c      archive complete. Clear archive buffers
c=======================================================================
c 
        elseif(msgtag.eq.MSG_ARC_CLR)then 
          call pvmfunpack( INTEGER_PVM, ittarch, 1, 1, info )
          if(ittarch.ne.ittsav) then
            write(outstr,901) ittarch, ittsav
            call s_send(MSG_OUTSTR)
 901  format(' Slave . Archive clear timestep mis-match',2i10)
          endif
c unset archive flag (in case next master request arrives early)
          archts = .false.
          larchp = .false.
          larchq = .false. 
          if(idebug.gt.3)then 
            write(ioslave,*)' Slave - message MSG_ARC_CLR received' 
          endif
c
c=======================================================================
c  9.  Receive message type MSG_TIMEVAR restart time variables
c=======================================================================
c
      elseif(msgtag.eq.MSG_TIMEVAR) then
        call pvmfunpack( INTEGER_PVM, itt   , 1, 1, info ) 
        call pvmfunpack( REAL_PVM  , totsec, 1, 1, info )
        call pvmfunpack( REAL_PVM  , totday, 1, 1, info )
        call pvmfunpack( REAL_PVM  , years , 1, 1, info )
        if(idebug.gt.2)then 
          write(ioslave,*) 
     &    ' Slave . Message MSG_TIMEVAR received.' ,itt
        endif
c
c=======================================================================
c 10.  Receive message type MSG_RESTART containing restart data
c=======================================================================
c
        elseif(msgtag.eq.MSG_RESTART)then
          call pvmfunpack( INTEGER_PVM,  idvar, 1, 1, info )
          call pvmfunpack( INTEGER_PVM, jarchl, 1, 1, info )
          call pvmfunpack( INTEGER_PVM, jarchu, 1, 1, info )
c
c convert master indices to slave indices
c
         jarchl = jarchl - jswm1
         jarchu = jarchu - jswm1
c
         if(idvar.le.NUM2D) then 
           do 60 j=jarchl,jarchu
           do 60 i=1,it
             call pvmfunpack(REAL_PVM, hsave(i,j,idvar),  1, 1, info) 
60         continue 
         elseif(vartgrd(idvar)) then 
           do 70 j=jarchl,jarchu
           do 70 i=1,it
             if(kmt(i,j).ne.0) call pvmfunpack(REAL_PVM, 
     &                     tsave(1,i,j,idvar-NUM2D), kmt(i,j), 1, info)
70         continue
         else
           do 80 j=jarchl,jarchu
           do 80 i=1,it
             if(kmu(i,j).ne.0) call pvmfunpack(REAL_PVM, 
     &                     tsave(1,i,j,idvar-NUM2D), kmu(i,j), 1, info)
80         continue
         endif 
c
c send a position message to indicate receipt of the data
c
         call sendposn(2)
c
          if(idebug.gt.3)then
            write(ioslave,*)' Slave - message MSG_RESTART received'
          endif
c
c=======================================================================
c 11.  Receive message type MSG_RQ_SNAP
c      timestep 'tsi' information request from master
c=======================================================================
c
        elseif(msgtag.eq.MSG_RQ_SNAP)then
          lsnapq = .true.
          call pvmfunpack( INTEGER_PVM, ittsnap2, 1, 1, info )
          call pvmfunpack( INTEGER_PVM, idsnp  , 1, 1, info )
c
          if(idebug.gt.3)then
            write(outstr,902) idsnp
            call s_send(MSG_OUTSTR)
            write(ioslave,*)' Slave - message MSG_RQ_SNAP received'
 902  format(' Slave . Message MSG_RQ_SNAP received, idsnp = ',i10)
          endif
c
c=======================================================================
c 12.  Receive message type MSG_BOUND1 
c      3-d boundary data from other slaves
c=======================================================================
c
      elseif (msgtag .eq. MSG_BOUND1) then
#ifdef pvm_buffer
        call pvmfunpack (REAL_PVM, bufout(1), 1, 1, info)
        lbuff = nint(bufout(1))
        if(lbuff.gt.LENBUF)then
           write(ioslave,*) ' error in s_when receiving MSG_BOUND2'
           write(ioslave,*) ' buffer overflow'
           call msgbye
        endif
        call pvmfunpack (REAL_PVM, bufout(2), lbuff-1, 1, info)
        mes   = nint(bufout(2))
        itts  = nint(bufout(3))
        ittps = nint(bufout(4))
        npts  = nint(bufout(5))
        lbuff = 5
#else
       call pvmfunpack( INTEGER_PVM, mes, 1, 1, info)
       call pvmfunpack( INTEGER_PVM, itts, 1, 1, info)
       call pvmfunpack( INTEGER_PVM, ittps, 1, 1, info)
       call pvmfunpack( INTEGER_PVM, npts, 1, 1, info)
#endif
c 
c       check information from this slave is expected 
c 
        do 85 na = 1,nactive 
         if(nactme(na).eq.mes) then
          if(bclnflg(na)) then 
           write(ioslave,*) 'Error duplicate or unexpected MSG_BOUND1' 
           write(ioslave,*) 
     &               'Receiving slave: ',me,' Sending slave: ',mes
           call s_send(MSG_BYE) 
          else 
           bclnflg(na) = .true. 
           goto 86 
          endif  
         endif  
   85   continue 
   86   continue
c
c       check timestep
c
        if(itts .ne. itt3) then
          write(ioslave,*) ' Error: Message MSG_BOUND1 from slave ',
     &                       mes,'has wrong timeflag'
          write(ioslave,*) ' itt from slave is ',itts, 
     &                     ' should be ',itt3
         call s_send(MSG_BYE)
        endif
       do 90 l = 1, npts
#ifdef pvm_buffer
          ii = nint(bufout(lbuff+1))
          jj = nint(bufout(lbuff+2))
          lbuff = lbuff+2
#else
        call pvmfunpack (INTEGER_PVM, ii, 1, 1, info)
        call pvmfunpack (INTEGER_PVM, jj, 1, 1, info)
#endif
        i = ii - iswm1
        if(i.le.0) i = i + IMT_M
        j = jj - jswm1
#ifdef pvm_buffer
          do 51 k=1, kmt(i,j)
            lbuff = lbuff+1
            t(k,i,j,1,np3) = bufout(lbuff)
  51      continue
          do 52 k=1, kmt(i,j)
            lbuff = lbuff+1
            t(k,i,j,2,np3) = bufout(lbuff)
  52      continue
          do 53 k=1, kmu(i,j)
            lbuff = lbuff+1
           u(k,i,j,np3) =  bufout(lbuff)
  53      continue
          do 54 k=1, kmu(i,j)
            lbuff = lbuff+1
            v(k,i,j,np3) = bufout(lbuff)
  54      continue
#else
        call pvmfunpack (REAL_PVM, t(1,i,j,1,np3), kmt(i,j), 1, info)
        call pvmfunpack (REAL_PVM, t(1,i,j,2,np3), kmt(i,j), 1, info)
        call pvmfunpack (REAL_PVM, u(1,i,j,np3), kmu(i,j), 1, info)
        call pvmfunpack (REAL_PVM, v(1,i,j,np3), kmu(i,j), 1, info)
#endif
#ifndef SYNC_3
c
c  if model is cyclic and point is in the end rows, copy it across
c
          if(icycl.eq.1)then
            if(i.eq.2)then
              do 541 k = 1,kmt(i,j)
                t(k,it,j,1,np3) = t(k,2,j,1,np3)
                t(k,it,j,2,np3) = t(k,2,j,2,np3)
 541          continue
               do 542 k = 1,kmt(i,j)
                u(k,it,j,np3) = u(k,2,j,np3)
                v(k,it,j,np3) = v(k,2,j,np3)
 542          continue
            elseif(i.eq.itm1)then
              do 543 k = 1,kmu(i,j)
                t(k,1,j,1,np3) = t(k,itm1,j,1,np3)
                t(k,1,j,2,np3) = t(k,itm1,j,2,np3)
 543          continue
              do 544 k = 1,kmu(i,j)
                u(k,1,j,np3) = u(k,itm1,j,np3)
                v(k,1,j,np3) = v(k,itm1,j,np3)
 544          continue
            endif
          endif
#endif
c
c       check data not already received and set flag
c
        if (iota(i,j,1)) then
          write (ioslave,*) 'Error: Message MSG_BOUND1 from slave',
     &    mes,' is overwriting previously received data point ', i,j
c         call s_send(MSG_BYE)
        else
           iota(i,j,1) = .true.
         nfalse(1) = nfalse(1) - 1
        endif
   90   continue
c
c check if all active neighbours have sent data
c
        ibcln = nactive
        do 94 na = 1,nactive
          if(bclnflg(na)) ibcln = ibcln - 1
   94   continue 
c
        if(idebug.gt.4)then
          write(ioslave,*)' Slave - message MSG_BOUND1 received.',
     &                    ' iposn=',iposn
        endif
        return
c
c=======================================================================
c 13.  Receive message type MSG_BOUND2 
c      2-d boundary data from other slaves
c=======================================================================
c
      elseif (msgtag .eq. MSG_BOUND2) then
#ifdef pvm_buffer
        call pvmfunpack (REAL_PVM, bufout(1), 1, 1, info)
        lbuff = nint(bufout(1))
        if(lbuff.gt.LENBUF)then
          write(ioslave,*) ' error in s_when receiving MSG_BOUND2'
          write(iosalve,*) ' buffer overflow'
          call msgbye
        endif
        call pvmfunpack (REAL_PVM, bufout(2), lbuff-1, 1, info)
        mes   = nint(bufout(2))
        itts  = nint(bufout(3))
        itbts = nint(bufout(4))
        itbtps= nint(bufout(5))
        npts  = nint(bufout(6))
        lbuff = 6
#else
        call pvmfunpack( INTEGER_PVM,    mes, 1, 1, info)
        call pvmfunpack( INTEGER_PVM,   itts, 1, 1, info)
        call pvmfunpack( INTEGER_PVM,  itbts, 1, 1, info)
        call pvmfunpack( INTEGER_PVM, itbtps, 1, 1, info)
        call pvmfunpack( INTEGER_PVM,   npts, 1, 1, info)
#endif
c
c       check information from this slave is expected
c
        do 95 na = 1,nactive
          if(nactme(na).eq.mes) then
            if(btrpflg(na)) then
              write(ioslave,*) 'Error duplicate or unexpected',
     &                         ' MSG_BOUND2'
              write(ioslave,*) 
     &               'Receiving slave: ',me,' Sending slave: ',mes
              call s_send(MSG_BYE)
            else
              btrpflg(na) = .true.
              goto 96
            endif
          endif
   95   continue
   96   continue
c
c       check timesteps
c
        if (itts .ne. itt2) then
          write(ioslave,*) 'Error: Message MSG_BOUND2 from slave ',mes,
     &                     'has wrong timeflag'
          write(ioslave,*) 'itt from slave is ',itts, 'should be ',
     &                      itt2
          call s_send(MSG_BYE)
        elseif (itbts .ne. itbt2) then
          write(ioslave,*) 'Error: Message MSG_BOUND2 from slave ',mes,
     &                     'has wrong free surface timeflag'
          write(ioslave,*) 'itbt from slave is ',itbts, 'should be ',
     &                      itbt2
          call s_send(MSG_BYE)
        elseif (itbtps .ne. itbtp2) then
          write(ioslave,*) 'Error: Message MSG_BOUND2 from slave ',mes,
     &                     'has wrong free surface pass number'
          write(ioslave,*) 'itbtp from slave is ',itbtps, 'should be ',
     &                      itbtp2
          call s_send(MSG_BYE)
        endif
c
        do 100 l = 1, npts
#ifdef pvm_buffer
          ii = nint(bufout(lbuff+1))
          jj = nint(bufout(lbuff+2))
#else
          call pvmfunpack (INTEGER_PVM, ii, 1, 1, info)
          call pvmfunpack (INTEGER_PVM, jj, 1, 1, info)
#endif
          i = ii - iswm1
          if(i.le.0) i = i + IMT_M
          j = jj - jswm1
#ifdef pvm_buffer
          h0(i,j,np02) = bufout(lbuff+3)
          u0(i,j,np02) = bufout(lbuff+4)
          v0(i,j,np02) = bufout(lbuff+5)
          lbuff = lbuff+5
#else
          call pvmfunpack (REAL_PVM, h0(i,j,np02), 1, 1, info)
          call pvmfunpack (REAL_PVM, u0(i,j,np02), 1, 1, info)
          call pvmfunpack (REAL_PVM, v0(i,j,np02), 1, 1, info)
#endif
#ifndef SYNC_2
c
c  if model is cyclic and point is in the end rows, copy it across
c
          if(icycl.eq.1)then
            if(i.eq.2)then
              h0(it,j,np02) = h0(2,j,np02)
              u0(it,j,np02) = u0(2,j,np02)
              v0(it,j,np02) = v0(2,j,np02)
            elseif(i.eq.itm1)then
              h0(1,j,np02) = h0(itm1,j,np02)
              u0(1,j,np02) = u0(itm1,j,np02)
              v0(1,j,np02) = v0(itm1,j,np02)
            endif
          endif
#endif
c
c       check data not already received and set flag
c
          if (iota(i,j,2)) then
            write (ioslave,'(a,i4,a,i4,a,2i4,a,2i4,a)')
     &     'Error: Message MSG_BOUND2 from slave ',mes,
     &     'is overwriting previously received data point in slave ', 
     &      me,' at ',i,j,'= Master: (',ii,jj,')'
c           call s_send(MSG_BYE)
          else
            iota(i,j,2) = .true.
            nfalse(2) = nfalse(2) - 1
          endif
  100   continue
c
c check if all active neighbours have sent data
c
        ibtrp = nactive
        do 105 na = 1,nactive
          if(btrpflg(na)) ibtrp = ibtrp - 1
  105   continue
c
        if(idebug.gt.4)then
          write(ioslave,*)' Slave - message MSG_BOUND2 received.',
     &                    ' iposn=',iposn
        endif
        return
c
c=======================================================================
c 14.  Receive message type MSG_CLR_BD1
c      Active neighbour ready to receive more baroclinic data
c=======================================================================
c
      elseif (msgtag .eq. MSG_CLR_BD1) then
        call pvmfunpack( INTEGER_PVM, mes, 1, 1, info)
c
        iclnclr = iclnclr + 1
        if(iclnclr.eq.nactive) lbcln = .true.
        if(idebug.gt.4)then
          write(ioslave,*)' Slave - message MSG_CLR_BD1 received.',
     &                    ' from ',mes
        endif
c
c
c=======================================================================
c 15.  Receive message type MSG_CLR_BD2
c      Active neighbour ready to receive more barotropic data
c=======================================================================
c
      elseif (msgtag .eq. MSG_CLR_BD2) then
        call pvmfunpack( INTEGER_PVM, mes, 1, 1, info)
c 
        itrpclr = itrpclr + 1
        if(itrpclr.eq.nactive) lbtrp = .true.
        if(idebug.gt.4)then
          write(ioslave,*)' Slave - message MSG_CLR_BD2 received.',
     &                    ' from ',mes
        endif
c 
c
c=======================================================================
c 16.  error
c=======================================================================
c
        else
          write(outstr,903)msgtag
          write(stderr,*) outstr
          call s_send(MSG_OUTSTR)
 903  format( ' Slave . Message type: ',i5,' not yet implemented')
      endif
c
c  if msgtype < 0 then jump back to process any remaining messages
c                 else return to calling subroutine
c
        if(msgtype.lt.0)goto 10
        return
      endif
c
c=======================================================================
c  PART 2.
c  send off any tsi, snapshot or archive messages now ready
c  to go
c=======================================================================
c
c  1.  send timestep information if ready
c
      if(ltsip .and. ltsiw)then
        call s_send(MSG_RY_TSI)
        ltsiw = .false.
        ltsip = .false.
      endif
c
c  2.  send partial archive data if ready
c
      if(larchp .and. larchq) then
        call s_send(MSG_RY_ARC)
        larchq=.false.
      endif
c
c
c  2.  send partial snapshot data if ready
c
      if(lsnapp .and. lsnapq) then
        call s_send(MSG_RY_SNAP)
        lsnapq=.false.
      endif
c
      return
      end
