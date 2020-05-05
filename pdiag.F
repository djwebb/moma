      subroutine pdiag
c
c @(#) SCCS module: pdiag.F  version: 1.3
c     Creation date: 03/17/95
c
c==================================================================
c     write diagnostics
c==================================================================
c
#include "def_slave.h"
#include "param.h"
#include "coord.h"
#include "ctmngr.h"
#include "iounit.h"
#include "frees.h"
#include "timelv.h"
#include "slabs.h"
      dimension ptemp(IMT_S,KM)
      character*12 trname(nt)
c
c==================================================================
c based on MOM common "cprnts.h" - only used by this routine
c
      parameter (nlatpr=10)
c
      dimension prlat(nlatpr), prslon(nlatpr), prelon(nlatpr)
     &,         prsdpt(nlatpr), predpt(nlatpr)
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
      jj = indp (prlat(jlat), yt, JMT_S)
      if (prlat(jlat) .le. yt(JMT_S)) then
        if (jlat .eq. 1) write(ioslave,9000)
        is = indp (prslon(jlat), xt, IMT_S)
        ie = indp (prelon(jlat), xt, IMT_S)
        ks = indp (prsdpt(jlat), zt, KM)
        ke = indp (predpt(jlat), zt, KM)
        fx = 1.0e-2
        do 680 m=1,nt
        scl = c1
        if (m .eq. 2) scl=1.e-3
        do 710 k=1,KM
        do 710 i=1,IMT_S
 710    ptemp(i,k) = t(k,i,jj,m,nc)
        write (ioslave,9100) trname(m), itt
     &,    jj, yt(jj), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
        call matrix (ptemp(1,1), IMT_S, is, ie, ks, ke, scl)
680     continue
c
        do 720 k=1,KM
        do 720 i=1,IMT_S
 720    ptemp(i,k) = u(k,i,jj,nc)
        scl = c1
        write (ioslave,9100) 'u velocity', itt
     &,   jj, yt(jj), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
        call matrix (ptemp(1,1), IMT_S, is, ie, ks, ke, scl)
c
        do 730 k=1,KM
        do 730 i=1,IMT_S
 730    ptemp(i,k) = v(k,i,jj,nc)
        scl = c1
        write (ioslave,9100) 'v velocity', itt
     &,   jj, yt(jj), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
        call matrix (ptemp(1,1), IMT_S, is, ie, ks, ke, scl)
      endif
700   continue
c
c-------------------------------------------------------------------
c       print the free surface model solution
c-------------------------------------------------------------------
c
      scl = c100
      is = indp (slonxy, xt, IMT_S)
      ie = indp (elonxy, xt, IMT_S)
      js = indp (slatxy, yt, JMT_S)
      je = indp (elatxy, yt, JMT_S)
      write (ioslave,8000) ' surface height (m)'
     &, itt, xt(is), xt(ie), yt(js), yt(je), scl
      call matrix (h0(1,1,np0), IMT_S, is, ie, -js, -je, scl)
c
      scl = c1
      write (ioslave,8000) ' u velocity (cm/s)'
     &, itt, xt(is), xt(ie), yt(js), yt(je), scl
      call matrix (u0(1,1,np0), IMT_S, is, ie, -js, -je, scl)
c
      write (ioslave,8000) ' v velocity (cm/s)'
     &, itt, xt(is), xt(ie), yt(js), yt(je), scl
      call matrix (v0(1,1,np0), IMT_S, is, ie, -js, -je, scl)
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
