      subroutine init_kmt
c
c @(#) SCCS module: init_kmt.F  version: 1.8
c     Creation date: 08/23/95
c
c==================================================================
c
c Subroutine to set-up the kmt and kmu arrays on the master
c
c==================================================================
c
#include "def_master.h"
#include "param.h"
#include "levind.h"
#include "coord.h"
#include "grdvar.h"
#include "switch.h"
#include "mesdta.h"
c
      character*8 string8
      logical llog1,llog2
c
c----------------------------------------------------------------------
c  read the kmt field either from ocean.kmt (if 'init' is set) ...
c----------------------------------------------------------------------
c
      if(init)then
        open(56,file='ocean.kmt',status='old',iostat = iostat)
        if(iostat.ne.0)then
          write(stdout,11)iostat
          stop
        endif
        read(56,15)string8,imt,jmt
        if(imt.ne.IMT_M.or.jmt.ne.JMT_M.or.string8.ne.'kmt     ')then
          if(string8.ne.'kmt     ')then
            write(stdout,12)string8
          else
            write(stdout,13)imt,jmt,IMT_M,JMT_M
          endif
          stop
        endif
        do 20 j=jmt,1,-1
          read(56,16)(kmt(i,j),i=1,imt)
  20    continue
        close(56)
c
c-----------------------------------------------------------------------
c  ...  or from the archive file
c-----------------------------------------------------------------------
c
      else
       call archrd(0)
      endif
c
c------------------------------------------------------------------
c     initialise kmu array
c------------------------------------------------------------------
c
      do 200 i=1,IMT_M
        kmu(i,JMT_M) = 0
 200  continue
      do 210 j=1,JMTM1_M
#ifdef cyclic_master
      do 210 i=1,IMT_M
        ip1 = mod(i,IMT_M) + 1
#else
      do 210 i=2,IMTM1_M
        ip1 = i+1
#endif
        kmu(i,j) = min(kmt(i,j),kmt(ip1,j),kmt(i,j+1),kmt(ip1,j+1))
 210  continue
c
c------------------------------------------------------------------
c     compute area and volume of ocean ("t,s" grid boxes)
c     and set cyclic boundary conditions
c------------------------------------------------------------------
c
      area   = c0
      volume = c0
      do 300 j=1,JMT_M
      do 300 i=1,IMT_M
        if (kmt(i,j) .gt. 0) then
          area   = area + cst(j)*dx*dy
          volume = volume + cst(j)*dx*dy*zw(kmt(i,j))
        endif
 300  continue
      write(stdout,301) area,volume
c
c------------------------------------------------------------------
c       print map of "kmt" levels.
c------------------------------------------------------------------
c
      if(IMT_M*JMT_M.le.8000)then
        write (stdout,401)
        do 410 ibk=1,IMT_M,40
          isp  = ibk  
          iept = ibk + 40 - 1
          if(iept.gt.IMT_M) iept=IMT_M
          write (stdout,'(/, 4x, 40i3)') (ii, ii=isp,iept)
          do 400 jrev=1,JMT_M
            j=JMT_M-jrev+1
            write (stdout,'(1x,i3, 40i3)')j,(kmt(i,j),i=isp,iept)
 400      continue    
 410    continue 
      else
        write(stdout,402)
      endif
c
c------------------------------------------------------------------
c    check kmt for consistency
c------------------------------------------------------------------
c
      llog1 = .false.
      llog2 = .false.
      do 450 i=1,IMT_M
        if(kmt(i,1).ne.0.or.kmt(i,JMT_M).ne.0)llog1=.true.
 450  continue
      if(llog1)then
        write(stdout,*)' Error: kmt array has non-zero values',
     &                 ' on top or bottom rows.'
        write(stdout,*)' Program stopping ...'
        stop
      endif
      do 460 j=2,JMTM1_M
        if(kmt(1,j).ne.0)llog1=.true.
        if(kmt(IMT_M,j).ne.0)llog2=.true.
 460  continue
#ifdef cyclic_master
      if(.not.llog1 .or. .not.llog2)then
        write(stdout,*)' Warning: flag "cyclic_master" set but ',
     &  'kmt array has only zero values in first or last column.'
        write(stdout,461)(kmt(1,j),j=2,JMTM1_M)
        write(stdout,461)(kmt(IMT_M,j),j=2,JMTM1_M)
 461  format(1x,'kmt: ',20i3)
         stop
      endif
#else
      if(llog1 .or. llog2)then
        write(stdout,*)' Error: flag "cyclic_master" not set but ',
     &  'kmt array has non-zero values in first or last column.'
        write(stdout,*)' Program stopping ...'
        stop
      endif
#endif
c
c------------------------------------------------------------------
c    set the sea-point counters for each block to be processed by
c    the master during archiving
c------------------------------------------------------------------
c
      nseat = 0
      nseau = 0
      do 510 n = 1,NZONE_M
        jlow = 1+(n-1)*JSUB_M
        jupp = min(JMT_M,n*JSUB_M)
        nsubst(n) = 0
        nsubsu(n) = 0
        do 500 j=jlow,jupp
        do 500 i=1,IMT_M
          if(kmt(i,j).ne.0) nsubst(n) = nsubst(n) + 1
          if(kmu(i,j).ne.0) nsubsu(n) = nsubsu(n) + 1
 500    continue
        nseat = nseat + nsubst(n)
        nseau = nseau + nsubsu(n)
 510  continue
      return
c
  11  format(1x,'Error opening file "ocean.kmt".  iostat = ',i4,/,
     &      ' Program stopping')
  12  format(1x,'Error opening file "ocean.kmt".  File type faulty',/,
     &      ' File type = ',a8,/,' Program stopping')
  13  format(1x,'Error opening file "ocean.kmt".  File dimensions ',
     &      'incompatable with current model.',/,
     &      ' imt, jmt from file =',2i5,/,' imt, jmt of model  =',2i5,/,
     &      ' Program stopping')
  15  format(1x,a8,6i10)
  16  format(1x,30i3)
 301  format(//,' Regional & Global ocean statistics (Master):'
     &,/,' the total ocean surface area (t grid)=',1pe15.8,'cm**2'
     &,/,' the total ocean volume (t grid)      =',1pe15.8,'cm**3')
 401  format(/t50,'number of levels on "t,s" grid')
 402  format(/t50,'number of levels on "t,s" grid',
     &             ' - array too large for listing')
      end
