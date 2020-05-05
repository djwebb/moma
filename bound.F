      subroutine bound
c
c @(#) SCCS module: bound.F  version: 1.8
c     Creation date: 08/23/95
c
c===================================================================
c
c    calculate boundary points for slave
c
c===================================================================
c
#include "def_slave.h"
#include "param.h"
#include "scalar.h"
#include "coord.h"
#include "iounit.h"
#include "levind.h"
#include "mesdta.h"
      logical log1
c
      do 20 i = 1, IBOUND
       ia(i) = 0
       ja(i) = 0
       do 10 k = 1,3
        isd(i,k) = 0
   10  continue
   20 continue
c
c     calculate interior points
c
      m = 0
      do 30 j = 2, jtm1
      do 30 i = 2, itm1
       if (kpn(i,j) .eq. me) then
        if ((kpn(i-1,j-1) .eq. me .or. kpn(i-1,j-1) .eq. 0) .and.
     &      (kpn(i-1,j  ) .eq. me .or. kpn(i-1,j  ) .eq. 0) .and.
     &      (kpn(i-1,j+1) .eq. me .or. kpn(i-1,j+1) .eq. 0) .and.
     &      (kpn(i  ,j+1) .eq. me .or. kpn(i  ,j+1) .eq. 0) .and.
     &      (kpn(i+1,j+1) .eq. me .or. kpn(i+1,j+1) .eq. 0) .and.
     &      (kpn(i+1,j  ) .eq. me .or. kpn(i+1,j  ) .eq. 0) .and.
     &      (kpn(i+1,j-1) .eq. me .or. kpn(i+1,j-1) .eq. 0) .and.
     &      (kpn(i  ,j-1) .eq. me .or. kpn(i  ,j-1) .eq. 0)) then
         m = m+1
         ia(m) = i
         ja(m) = j
        endif
       endif
   30 continue
      ninner = m
      write (ioslave,*) 'number of interior points: ', ninner
c
c     calculate inner halo points
c
      m = 0
      do 50 j = 2, jtm1
      do 50 i = 2, itm1
       log1 = .false.
       if (kpn(i,j) .eq. me) then
        do 40 jj = j-1, j+1
        do 40 ii = i-1, i+1
         if (kpn(ii,jj) .ne. me .and. kpn(ii,jj) .ne .0) log1 = .true.
   40    continue
         if (log1) then
          m = m + 1
          ia(ninner+m) = i
          ja(ninner+m) = j
         endif
       endif
   50 continue
      nhalo = m
      ncore = nhalo + ninner
      write (ioslave,*) 'number of inner halo points: ', nhalo
      write (ioslave,*) 'total number of sea points: ', ncore
c
      m = 0
      do 80 l = 1, nhalo
       k = ninner + l
       ik = ia(k)
       jk = ja(k)
       log1 = .false.
       do 70 j = jk-1,jk+1
       do 70 i = ik-1,ik+1
        if (kpn(i,j) .ne. me .and. kmt(i,j) .ne. 0) then
         log1 = .true.
         if (icycl .eq. 1) then
          if ((i .eq. 1 .or. i .eq. it) .and. kpn(1,j) .eq. kpn(it,j))
     &     log1 = .false.
         endif
         if (log1) then
          if (m .eq. 0) then
           m = m+1
           isd(m,1) = ik
           isd(m,2) = jk
           isd(m,3) = kpn(i,j)
          else
           do 60 jj = 1, m
            if (isd(jj,1) .eq. ik .and. isd(jj,2) .eq. jk .and.
     &          isd(jj,3) .eq. kpn(i,j)) log1 = .false.
   60      continue
           if (log1) then
            m=m+1
            isd(m,1) = ik
            isd(m,2) = jk
            isd(m,3) = kpn(i,j)
            if (m .gt. IBOUND-1) then
             print *, 'error - too many boundary points'
             stop
            endif
           endif
          endif
         endif
        endif
   70  continue
   80 continue
      nbound = m
      write (ioslave,*) 'number of boundary points: ', nbound
c
c     sort according to neighbouring processor numbers
c
      nmax = 0
      do 90 i = 1, nbound
       if (isd(i,3) .gt. nmax) nmax = isd(i,3)
   90 continue
      write (ioslave, *) 'nmax: ',nmax
c
      m = 0
      do 110 n = nmax,1,-1
      do 110 i = 1, nbound
       if (isd(i,3) .eq. n) then
        m = m+1
        do 100 k = 1,3
         itemp = isd(i,k)
         isd(i,k) = isd(m,k)
         isd(m,k) = itemp
  100   continue
       endif
  110 continue
c
      do 120 n = 1, nproc
       nsend(n) = 0
  120 continue
      do 130 k = 1,nbound
       l = isd(k,3)
       nsend(l) = nsend(l)+1
  130 continue
      write (ioslave,*) 'nsend: ',(i, nsend(i),i=1,nmax)
c
c store the identifiers of the active neighbours (i.e.
c those neighbours which receive information from me
c
      nactive=0
      do 135 n=1,nproc
       if(nsend(n).gt.0) then
        nactive=nactive+1
        nactme(nactive)=n
        nactid(nactive)=itida(n)
       endif
  135 continue
c
c     calculate outer halo points
c
      m = 0
      do 150 j = 1, jt
      do 150 i = 1, it
       if (kpn(i,j) .ne. me .and. kmt(i,j) .ne. 0) then
        log1 = .false.
        ip1 = min (i+1,IMT_S)
        im1 = max (1,i-1)
        jp1 = min (j+1,JMT_S)
        jm1 = max (j-1,1)
        do 140 l = jm1,jp1
        do 140 k = im1,ip1
         if (kpn(k,l) .eq. me) log1 = .true.
  140   continue
        if (icycl .eq. 1) then
         if (i .eq. 1 .or. i .eq. it)
     &     log1 = .false.
        endif
        if (log1) then
         m = m+1
         iout(m) = i
         jout(m) = j
        endif
       endif
  150 continue
      nouter = m
      write (ioslave,*) 'nouter: ', nouter  
      return
      end
