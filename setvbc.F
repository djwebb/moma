      subroutine setvbc (ic, jc)
c
c @(#) SCCS module: setvbc.F  version: 1.9
c     Creation date: 12/18/95
c
c==================================================================
c
c     compute vertical boundary conditions for clinic & tracer
c
c==================================================================
c
#include "def_slave.h"
#include "param.h"
#include "cvbc.h"
#include "coord.h"
#include "ctmngr.h"
#include "grdvar.h"
#include "scalar.h"
#include "slabs.h"
#include "timelv.h"
#include "levind.h"
#include "iounit.h"
c
c----------------------------------------------------------------
c  1.  interpolate wind stress data
c   use linear interpolation to produce the estimated surface
c   boundary condition values for wind stress components at
c   velocity point (ic,jc)
c----------------------------------------------------------------
c
      xlon = xu(ic)
      ylat = yu(jc)
c
c----------------------------------------------------------------
c   first find wind stress indices south-west of velocity point
c----------------------------------------------------------------
c
      xx = xlon - stlon_w
      yy = ylat - stlat_w
      ii = xx/dx_w
      jj = yy/dy_w
c
      if(lcycw)then
        if(xx.lt.0)ii = ii - 1
        ii = mod(ii,i_w)+1
        ip = ii+1
        if(ip.gt.i_w)ip = 1
      else
        ii = ii + 1
        if(ii.eq.i_w)ii=i_w-1
        ip = ii + 1
      endif
c
      jj = jj + 1
      if(jj.eq.j_w)jj=j_w-1
      jp = jj+1
      if(xx.lt.0.0.or.ii.gt.i_w .or.
     &   yy.lt.0.0.or.jj.gt.j_w )then
          write(ioslave,901)ic,jc,ii,jj,i_w,j_w,xx,yy,xlon,ylat,
     &                      stlon_w,stlat_w,dx_w,dy_w,lcycw
 901  format(' Routine setvbc: error while finding wind indices',
     &       /,'  ic, jc, ii, jj, i_w, j_w, xx, yy: ',6i5,2f10.4,
     &       /,'  xlon, ylat, stlon_w, stlat_w,',
     &         ' dx_w, dy_w, lcycw:',6f10.4,L)
         call s_abort
      endif
c
c----------------------------------------------------------------
c   interpolate in space and time
c----------------------------------------------------------------
c
      nmm = mod(mlast, 3) + 1
      nmp = mod(mnext, 3) + 1
c
c   calculate spatial weighting coefficients
c
      xwa = (xlon - (stlon_w + (ii-1)*dx_w))/dx_w
      xwb = 1.0 - xwa
      ywa = (ylat - (stlat_w + (jj-1)*dy_w))/dy_w
      ywb = 1.0 - ywa
c
c   interpolate x-component in space and time
c
      wupm = xwa*wstx(ip,jp,nmm) + xwb*wstx(ii,jp,nmm)
      wlom = xwa*wstx(ip,jj,nmm) + xwb*wstx(ii,jj,nmm)
      wupp = xwa*wstx(ip,jp,nmp) + xwb*wstx(ii,jp,nmp)
      wlop = xwa*wstx(ip,jj,nmp) + xwb*wstx(ii,jj,nmp)
c
      smf(1) = twa*(ywa * wupp + ywb * wlop) +
     &         twb*(ywa * wupm + ywb * wlom)
c       
c   interpolate y-component in space and time
c
      wupm = xwa*wsty(ip,jp,nmm) + xwb*wsty(ii,jp,nmm)
      wlom = xwa*wsty(ip,jj,nmm) + xwb*wsty(ii,jj,nmm)
      wupp = xwa*wsty(ip,jp,nmp) + xwb*wsty(ii,jp,nmp)
      wlop = xwa*wsty(ip,jj,nmp) + xwb*wsty(ii,jj,nmp)
c       
      smf(2) = twa*(ywa * wupp + ywb * wlop) +
     &         twb*(ywa * wupm + ywb * wlom)
c
c----------------------------------------------------------------
c  2.  interpolate met data
c   use linear interpolation to produce the estimated surface
c   boundary condition values for met data at tracer point (ic,jc)
c----------------------------------------------------------------
c
      xlon = xt(ic)
      ylat = yt(jc)
c
c----------------------------------------------------------------
c   first find sst indices south-west of tracer point
c----------------------------------------------------------------
c
      xx = xlon - stlon_s
      yy = ylat - stlat_s
      ii = xx/dx_s
      jj = yy/dy_s
c
      if(lcycw)then
        if(xx.lt.0)ii = ii - 1
        ii = mod(ii,i_s)+1
        ip = ii+1
        if(ip.gt.i_s)ip = 1
      else
        ii = ii + 1
        if(ii.eq.i_s)ii=i_s-1
        ip = ii+1
      endif
c
      jj = jj + 1
      if(jj.eq.j_s)jj=j_s-1
      jp = jj + 1
c
      if(xx.lt.0.0.or.ii.gt.i_s .or.
     &   yy.lt.0.0.or.jj.gt.j_s )then
          write(ioslave,902)ic,jc,ii,jj,i_s,j_s,xx,yy,xlon,ylat,
     &                      stlon_w,stlat_w,dx_w,dy_w,lcycw
 902  format(' Routine setvbc: error while finding flux indices',
     &       /,'  ic, jc, ii, jj, i_s, j_s, xx, yy: ',6i5,2f10.4,
     &       /,'  xlon, ylat, stlon_w, stlat_w,',
     &         ' dx_w, dy_w, lcycw:',6f10.4,L)
         call s_abort
      endif
c
c----------------------------------------------------------------
c   interpolate in space and time
c----------------------------------------------------------------
c
c   calculate spatial weighting coefficients
c
      xwa = (xlon - (stlon_s + (ii-1)*dx_s))/dx_s
      xwb = 1.0 - xwa
      ywa = (ylat - (stlat_s + (jj-1)*dy_s))/dy_s
      ywb = 1.0 - ywa
c
c   interpolate in space and time
c
      wupm = xwa*ssta(ip,jp,nmm) + xwb*ssta(ii,jp,nmm)
      wlom = xwa*ssta(ip,jj,nmm) + xwb*ssta(ii,jj,nmm)
      wupp = xwa*ssta(ip,jp,nmp) + xwb*ssta(ii,jp,nmp)
      wlop = xwa*ssta(ip,jj,nmp) + xwb*ssta(ii,jj,nmp)      
      stf(1) = twa*(ywa * wupp + ywb * wlop) +
     &         twb*(ywa * wupm + ywb * wlom)
c
      wupm = xwa*sssa(ip,jp,nmm) + xwb*sssa(ii,jp,nmm)
      wlom = xwa*sssa(ip,jj,nmm) + xwb*sssa(ii,jj,nmm)
      wupp = xwa*sssa(ip,jp,nmp) + xwb*sssa(ii,jp,nmp)
      wlop = xwa*sssa(ip,jj,nmp) + xwb*sssa(ii,jj,nmp)      
      stf(2) = twa*(ywa * wupp + ywb * wlop) +
     &         twb*(ywa * wupm + ywb * wlom)
c
c  convert salinity to model units.
c
      stf(2) = (stf(2)-35.0)*0.001
c
c------------------------------------------------------------------
c     set no flux condition for tracer at bottom.
c     relax to "observed" values with 50 day timescale at top
c------------------------------------------------------------------
c
      gamma = secday/50.0
      do 100 m=1,nt
        stf(m) = gamma*(stf(m)-t(1,ic,jc,m,nm))*dz(1)
        btf(m) = c0
 100  continue
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
