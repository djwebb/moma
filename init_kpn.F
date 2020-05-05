      subroutine init_kpn
c
c @(#) SCCS module: init_kpn.F  version: 1.10
c      Creation date: 08/09/96
c
c==================================================================
c
c  routine to read in kpn field (and kmt field if initialising)
c  and to initialise the arrays describing the start and end points
c  of each processor's region
c
c==================================================================
c
#include "def_master.h"
#include "param.h"
#include "coord.h"
#include "levind.h"
#include "mesdta.h"
#include "switch.h"
      logical*4 found
      integer incol(IMT_M), inrow(JMT_M)
      character*8 string8
c
c----------------------------------------------------------------------
c   read the processor map kpn from 'ocean.kpn'
c----------------------------------------------------------------------
c
      open(56,file='ocean.kpn',status='old',iostat = iostat)
      if(iostat.ne.0)then
        write(stdout,31)iostat
        stop
      endif
      read(56,15)string8,imt,jmt,nproc
      if(imt.ne.IMT_M.or.jmt.ne.JMT_M.or.string8.ne.'kpn     ')then
        if(string8.ne.'kpn     ')then
          write(stdout,32)string8
        else
          write(stdout,33)imt,jmt,IMT_M,JMT_M
        endif
        stop
      endif
c
      do 30 j=JMT_M,1,-1
       read(56,16) (kpn(i,j),i=1,IMT_M)
  30  continue
      close (56)
  31  format(1x,'Error opening file "ocean.kpn".  iostat = ',i4,/,
     &      ' Program stopping')
  32  format(1x,'Error opening file "ocean.kpn".  File type faulty',/,
     &      ' File type = ',a8,/,' Program stopping')
  33  format(1x,'Error opening file "ocean.kpn".  File dimensions ',
     &      'incompatable with current model.',/,
     &      ' imt, jmt from file =',2i5,/,' imt, jmt of model  =',2i5,/,
     &      ' Program stopping')
  15  format(1x,a8,6i10)
  16  format(1x,30i3)
c
c-----------------------------------------------------------------------
c  scan kpn and kmt arrays to check they are consistant
c-----------------------------------------------------------------------
c
      ierr1 = 0
      ierr2 = 0
      ierr3 = 0
      ierr4 = 0
      ierr5 = 0
      nproc1 = -999
      do 100 j=1,JMT_M
      do 100 i=1,IMT_M
        if(kmt(i,j).lt.0.or.kmt(i,j).gt.km)then
          ierr1 = ierr1+1
          if(ierr1.lt.10)write(stdout,*) ' Subroutine init_kpn: ',
     &     'kmt out of range.  i, j, kmt = ',i,j,kmt(i,j)
        endif
        if(kpn(i,j).lt.0.or.kpn(i,j).gt.nproc)then
          ierr2 = ierr2+1
          if(ierr2.lt.10)write(stdout,*) ' Subroutine init_kpn: ',
     &     'kpn out of range.  i, j, kpn, nproc = ',i,j,kmt(i,j),nproc
        endif
        if(kmt(i,j).ge.0.and.kmt(i,j).le.km.and.
     &    kpn(i,j).ge.0.and.kpn(i,j).le.nproc)then
          if(kmt(i,j).eq.0.and.kpn(i,j).ne.0)then
            ierr3 = ierr3+1
            if(ierr3.lt.10)write(stdout,*) ' Routine init_kpn: ',
     &       'kpn non zero when kmt zero.  i, j, kmt, kpn =',
     &       i,j,kmt(i,j),kpn(i,j)
          endif
          if(kpn(i,j).eq.0.and.kmt(i,j).ne.0)then
            ierr4 = ierr4+1
            if(ierr4.lt.10)write(stdout,*) ' Routine init_kpn: ',
     &      'kmt non zero when kpn zero.  i, j, kmt, kpn =',
     &       i,j,kmt(i,j),kpn(i,j)
          endif
        endif
      nproc1 = max(nproc1,kpn(i,j))
 100  continue
      if(ierr1+ierr2+ierr3+ierr4.ne.0.or.nproc.ne.nproc1)then
        write(stdout,*)' Routine init_kpn:'
        if(ierr1.ne.0)write(stdout,*)
     &    '  Number of kmt points out of range = ',ierr1
        if(ierr2.ne.0)write(stdout,*)
     &    '  Number of kpn points out of range = ',ierr2
        if(ierr3.ne.0)write(stdout,*)
     &    '  Number of kmt points set where kpn is zero = ',ierr3
        if(ierr4.ne.0)write(stdout,*)
     &    '  Number of kpn points set where kmt is zero = ',ierr4
        write(stdout,*) ' Program stopping '
        if(nproc.ne.nproc1)then
          write(stdout,*)'  nproc in ocean.kpn file header ',
     &    'inconsistant with array values'
          write(stdout,*)'  nproc =',nproc
          write(stdout,*)'  maximum value of kpn(i,j) = ',nproc1
        endif
        write(stdout,*) ' Program stopping '
        stop
      endif
c
c-----------------------------------------------------------------------
c  find the range spanned by each processor and set cyclic flag
c-----------------------------------------------------------------------
c
      do 200 n=1,nproc
        it_s(n)=2*imt_m
        it_e(n)=0
        jt_s(n)=jmt_m+1
        jt_e(n)=0
        icycl_M(n)=0
 200  continue 
c
c-----------------------------------------------------------------------
c  start loop over processors      
c-----------------------------------------------------------------------
c
      do 500 n=1,nproc
        found = .false.
        do 300 i=1,IMT_M
        incol(i) = 0
        do 300 j=1,JMT_M
          if(kpn(i,j).eq.n) then
            incol(i) = 1
            found = .true.
          endif
 300    continue
        if(.not.found) then
          write(stdout,*) ' Routine init_kpn:  faulty processor map'
          write(stdout,*) ' No points found for processor number: ', n
          write(stdout,*) ' Program stopping'
          stop
        endif
c
c-----------------------------------------------------------------------
c check which rows contain this processor
c-----------------------------------------------------------------------
c
        do 400 j=1,JMT_M
        inrow(j) = 0 
        do 400 i=1,IMT_M   
          if(kpn(i,j).eq.n) inrow(j) = 1
 400    continue
c
c-----------------------------------------------------------------------
c    Find the minimum spanning length of all the 1's in incol
c-----------------------------------------------------------------------
c
#ifdef cyclic_master
        call minspan(incol,IMT_M,it_s(n),it_e(n),icycl_M(n),1)
#else
        call minspan(incol,IMT_M,it_s(n),it_e(n),icycl_M(n),0)
        if(icycl_M(n).lt.0) 
     &       write(stdout,*) ' Routine init_kpn: i-span warning given',
     &                       ' for processor ',n
#endif
c
c-----------------------------------------------------------------------
c    Find the minimum spanning length of all the 1's in inrow
c-----------------------------------------------------------------------
c
        call minspan(inrow,JMT_M,jt_s(n),jt_e(n),idummy,0)
        if(idummy.lt.0) 
     &       write(stdout,*) ' Routine init_kpn: j-span warning given',
     &                       ' for processor ',n
c
c-----------------------------------------------------------------------
c     swlo and swla contain the long-lat co-ordinates of the sw
c     corner of the the most south-westerly halo point.
c-----------------------------------------------------------------------
c
        swlo(n) = stlon + (it_s(n)-1)*dxdeg
        swla(n) = stlat + (jt_s(n)-1)*dydeg
c
c-----------------------------------------------------------------------
c     iswm1 and jswm1 are the master indices for the cell to the 
c     south-west of the most south-westerly halo point (even it if
c     doesn't exist). These values are used as offsets by the master 
c     and slaves to determine the local indices of data passed 
c     through messages. I.e.:
c
c     slave(i_s,j_s) --> master(mod(iswm1+i_s,imt_m),jswm1+j_s)
c
c     Thus all messages passed from and to slaves which refer to 
c     specific cells should use the master indexing.
c-----------------------------------------------------------------------
c
        iswm1_M(n) = it_s(n)-1
        if (icycl_M(n) .eq. 1) iswm1_M(n) = -1
        jswm1_M(n) = jt_s(n)-1
        it_l(n) = it_e(n) - it_s(n) + 1
        jt_l(n) = jt_e(n) - jt_s(n) + 1
 500  continue
c
c-----------------------------------------------------------------------
c Write out slave bounding regions 
c-----------------------------------------------------------------------
c
      write(stdout,551) 'slave','i-start','i-end','i-length',
     &                            'j-start','j-end','j-length'
      do 550 n=1,nproc
        write(stdout,552) n,it_s(n),it_e(n),it_l(n),
     &                      jt_s(n),jt_e(n),jt_l(n)
 550  continue
      write(stdout,*)
c
 551  format(1x,/,7a10)
 552  format(7i10)
      return
      end


      subroutine minspan(intarr, ilenarr, istart, iend, kcycl, iwrap)
c
c subroutine to find the minimum spanning length of all the 1's in a
c binary integer array.
c
c intarr  - integer array containing sequences of 1's and 0's
c ilenarr - length of the integer array
c istart  - the start index of the minimum span - 1, for outer halo.
c iend    - the end index of the minimum span + 1, for outer halo.
c kcycl   - (1,0) if the span (is, is not) cyclic
c iwrap   - (1,0) if intarr (should, should not) be treated cyclically.
c
      integer intarr(ilenarr)
c
c
       minspn = ilenarr + 3
       istart = 0
       iend = 0 
       kcycl = 0
c

      if(iwrap.ne.0) then
c
       do 10 kstart = 1,ilenarr
        if(intarr(kstart).eq.0) goto 10
        lastone = kstart
        do 20 kk = kstart+1, kstart + ilenarr - 1
         k = mod(kk-1,ilenarr) + 1
         if(intarr(k).eq.1) lastone = kk
20      continue
        thispan = lastone - kstart + 1
        if(thispan.lt.minspn) then
         kcycl = 0
         minspn=thispan
         istart = kstart - 1
         if(istart.eq.0) istart = ilenarr
         iend = istart + minspn + 1
         if(istart.eq.ilenarr.and.iend.eq.2*ilenarr + 1) kcycl = 1
        endif
10     continue
c
      else
c
       do 30 kstart = 1,ilenarr 
        if(intarr(kstart).ne.0) goto 35  
30     continue
c
35      lastone = kstart 
        do 40 kk = kstart+1, ilenarr  
         if(intarr(kk).eq.1) lastone = kk   
40      continue  
        thispan = lastone - kstart + 1 
        if(thispan.lt.minspn) then  
         minspn=thispan 
         istart = kstart - 1
         if(istart.eq.0) istart = ilenarr
         iend = istart + minspn + 1    
         if(istart.eq.ilenarr.and.iend.eq.2*ilenarr + 1) then
          write(6,*) ' Warning: cyclic dimensions derived for'
          write(6,*) '          explicit non-cyclic domain.  '
          write(6,*) '          Possible error in kpn field. '
          kcycl = -1
         endif
        endif
c
      endif
c
      return
      end
