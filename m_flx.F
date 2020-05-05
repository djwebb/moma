      subroutine initmet
# ifdef hcomments
c
c @(#) SCCS module: m_flx.F  version: 1.9
c     Creation date: 12/18/95
c
c=======================================================================
c  master subroutine to initialise surface flux logic
c=======================================================================
c
c  This routine initialises the variables used to specify the surface
c  flux data sets.
c
c  The example implemented uses the default constant met data from the
c  mom model:
c
c   "observed" temperature and salinity data are based on global,
c   annual mean zonally averaged values from the Levitus Atlas (1982).
c
c   "observed" windstress data are based on global, annual mean,
c   zonally averaged values from Hellerman and Rosenstein (1981).
c   some smoothing was done.
c
c   Units:  temperature - degrees Centigrade
c           salinity    - standard oceanographic salinity units
c           wind stress - dynes per sq cm
c   The salinity is transformed into model units in routine setvbc'.
c
c   References:
c     Hellerman, S, and M. Rosenstein, normal monthly wind stress
c         over the world ocean with error estimates, Journal of
c         Physical Oceanography, 13, 1093-1104,1983.
c     Levitus, S., Climatological atlas of the world ocean, NOAA
c         Prof. Paper 13, US Gov't printing Office, Washington, DC, 
c         1982.
c
# endif
#include "def_master.h"
#include "param.h"
#include "coord.h"
#include "cvbc.h"
#include "mesdta.h"
c
      i_w     =  80
      j_w     =  NOLATP
      stlon_w =   4.5
      stlat_w = -90.0
      dx_w    =   4.5
      dy_w    =   4.5
c
      i_s     =  80
      j_s     =  NOLAT
      stlon_s =   2.25
      stlat_s = -87.75
      dx_s    =   4.5
      dy_s    =   4.5
c
c  check if the met fields are cyclic
c
      lcycw = abs(360.0 - dx_w*NWX_M).lt.0.01*dx_w
      lcycs = abs(360.0 - dx_s*NSX_M).lt.0.01*dx_s
      lcycm = lcycw.and.lcycs
c
c  loop over the processors finding the range each requires
c
      do 100 n=1,nproc
c
c  find latitude range for winds
c
        y1 = swla(n)-stlat_w
        j1 = max(1,1+int(y1/dy_w))
c
        if(y1.lt.-0.0001*dy_w)then
          write(stdout,901)'j1',n
          write(stdout,902)j1,y1,swla(n),stlat_w,dy_w
 901  format(' routine initmet: error while calculating variable ',
     &       a,' for slave ',i5)
 902  format(' j1 = ',i5,' y1 = ',f7.2,' swla(n)= ',f7.2,
     &       ' stlat_w = ',f7.2,' dy_w = ',f7.2)
        endif
c
        y2 = y1+jt_l(n)*dydeg
        j2 = 1+int(y2/dy_w)
c
        if(j2.lt.j_w)then
          if((j2-1)*dy_w.lt.y2)j2 = j2+1
        elseif(j2.eq.j_w .and.
     &      (j2-1)*dy_w.lt.y2-0.0001*dy_w)then
            write(stdout,901)'j2',n
            write(stdout,903)j1,j2,j_w,jt_l(n)
            write(stdout,904)y2,y1,dydeg
  903  format(' j1 = ',i5,' j2 = ',i5,' j_w = ',i5,' jt_l(n) = ',i5)
  904  format(' y2 = ',f7.2,' y1 = ',f7.2,' dydeg = ',f7.2)
       endif
c
        if(j1.lt.1.or.j1.gt.j_w.or.j2.lt.1.or.j2.gt.j_w)then
          write(stdout,901)'j1 and j2',n
          write(stdout,905)j1,j2,j_w,jt_l(n)
          write(stdout,906)swla(n),stlat_w,dy_w,dydeg
  905  format(' j1 = ',i5,' j2 = ',i5,' j_w = ',i5,' jt_l(n) =',i5)
  906  format(' swla(n) = ',f7.2,' stlat_w = ',f7.2,' dy_w = ',f7.2,
     &        ' dydeg = ',f7.2)
        endif
c
c  store latitude data for this slave
c
        jw_s(n) = j1
        jw_e(n) = j2
        jw_l(n) = j2-j1+1
        swla_w(n) = stlat_w+(j1-1)*dy_w
c
c  find longitude range for winds
c
        x1 = swlo(n)-stlon_w
        if(lcycm)then
          x1 = mod(x1,360.0)
          if(x1.lt.0.0)x1 = x1 + 360.0
        endif
        i1 = 1+int(x1/dx_w)
c
        if(.not.lcycm.and.x1.lt.-0.0001*dx_w)then
          write(stdout,901)'i1',n
          write(stdout,907)i1,x1,swlo(n),stlon_w,dx_w
 907  format(' i1 = ',i5,' x1 = ',f7.2,' swlo(n)= ',f7.2,
     &       ' stlon_w = ',f7.2,' dx_w = ',f7.2)
        endif
c
        x2 = x1+it_l(n)*dxdeg
        i2 = 1+int(x2/dx_w)
c
        if(lcycm)then
          if((i2-1)*dx_w.lt.x2)i2 = i2+1
          if(i2-i1.ge.i_w)   i2 = i1 +i_w-1
        elseif(i2.lt.i_w)then
          if((i2-1)*dx_w.lt.x2)i2 = i2+1
        elseif(i2.eq.i_w .and.
     &        (i2-1)*dx_w.lt.x2-0.0001*dx_w)then
          write(stdout,901)'i2',n
          write(stdout,908)i1,i2,i_w,it_l(n)
          write(stdout,909)x1,x2,swlo(n),stlon_w,dx_w,dxdeg
  908  format(' i1 = ',i5,' i2 = ',i5,' i_w = ',i5,' it_l(n) = ',i5)
  909  format(' x1 = ',f7.2,' x2 = ',f7.2,' swlo(n) = ',f7.2,
     &        ' stlon_w = ',f7.2,' dx_w = ',f7.2,' dxdeg = ',f7.2)
        endif
c
        iw_s(n) = i1
        iw_e(n) = i2
        iw_l(n) = i2-i1+1
        swlo_w(n) = stlon_w+(i1-1)*dx_w
c
c  find latitude range for surface values
c
        y1 = swla(n)-stlat_s
        j1 = max(1,1+int(y1/dy_s))
c
        if(y1.lt.-0.0001*dy_s)then
          write(stdout,901)n
          write(stdout,912)j1,y1,swla(n),stlat_s,dy_s
 912  format(' j1 = ',i5,' y1 = ',f7.2,' swla(n)= ',f7.2,
     &       ' stlat_s = ',f7.2,' dy_s = ',f7.2)
        endif
c
        y2 = y1+jt_l(n)*dydeg
        j2 = 1+int(y2/dy_s)
c
        if(j2.lt.j_s)then
          if((j2-1)*dy_s.lt.y2)j2 = j2+1
        elseif(j2.eq.j_s .and.
     &      (j2-1)*dy_s.lt.y2-0.0001*dy_s)then
            write(stdout,901)'j2',n
            write(stdout,913)j1,j2,j_s,jt_l(n)
            write(stdout,914)y2,y1,dydeg
  913  format(' j1 = ',i5,' j2 = ',i5,' j_s = ',i5,' jt_l(n) = ',i5)
  914  format(' y2 = ',f7.2,' y1 = ',f7.2,' dydeg = ',f7.2)
        endif
c
        if(j1.lt.1.or.j1.gt.j_s.or.j2.lt.1.or.j2.gt.j_s)then
          write(stdout,901)'j1 and j2',n
          write(stdout,915)j1,j2,j_s,jt_l(n)
          write(stdout,916)swla(n),stlat_s,dy_s,dydeg
  915  format(' j1 = ',i5,' j2 = ',i5,' j_s = ',i5,' jt_l(n) =',i5)
  916  format(' swla(n) = ',f7.2,' stlat_s = ',f7.2,' dy_s = ',f7.2,
     &        ' dydeg = ',f7.2)
        endif
c
c  store latitude data for this slave
c
        js_s(n) = j1
        js_e(n) = j2
        js_l(n) = j2-j1+1
        swla_s(n) = stlat_s+(j1-1)*dy_s
c
c  find longitude range for surface values
c
        x1 = swlo(n)-stlon_s
        if(lcycm)then
          x1 = mod(x1,360.0)
          if(x1.lt.0.0)x1 = x1 + 360.0
        endif
        i1 = 1+int(x1/dx_s)
c
        if(.not.lcycm.and.x1.lt.-0.0001*dx_s)then
          write(stdout,901)'x1',n
          write(stdout,917)i1,x1,swlo(n),stlon_s,dx_s
 917  format(' i1 = ',i5,' x1 = ',f7.2,' swlo(n)= ',f7.2,
     &       ' stlon_s = ',f7.2,' dx_s = ',f7.2)
        endif
c
        x2 = x1+it_l(n)*dxdeg
        i2 = 1+int(x2/dx_s)
c
        if(lcycm)then
          if((i2-1)*dx_s.lt.x2)i2 = i2+1
          if(i2-i1.ge.i_s)   i2 = i1 +i_s-1
        elseif(i2.lt.i_s)then
          if((i2-1)*dx_s.lt.x2)i2 = i2+1
        elseif(i2.eq.i_s .and.
     &        (i2-1)*dx_s.lt.x2-0.0001*dx_s)then
          write(stdout,901)'i2',n
          write(stdout,918)i1,i2,i_s,it_l(n)
          write(stdout,919)x1,x2,swlo(n),stlon_s,dx_s,dxdeg
  918  format(' i1 = ',i5,' i2 = ',i5,' i_w = ',i5,' it_l(n) = ',i5)
  919  format(' x1 = ',f7.2,' x2 = ',f7.2,' swlo(n) = ',f7.2,
     &        ' stlon_s = ',f7.2,' dx_s = ',f7.2,' dydeg = ',f7.2)
        endif
c
        is_s(n) = i1
        is_e(n) = i2
        is_l(n) = i2-i1+1
        swlo_s(n) = stlon_s+(i1-1)*dx_s
c
 100  continue
      return
c
      end


      subroutine metrd(month)
#ifdef hcomments
c
c @(#) SCCS module: m_flx.F  version: 1.9
c    Creation date: 12/18/95
c
c=======================================================================
c subroutine to read in the meteorological data 
c=======================================================================
c
c  This routine inputs the data defining the surface fluxes for month
c  'month' and stores it in arrays ready for transfer to the slave 
c  processes.  
c
c  The example implemented uses the default constant mom data, here
c  stored in arrays setup by routine initmet.  In normal use the data
c  would be read from disk as required.
c
c  The present model is set up to use data centered on the middle of
c  each month.  It should be straightforward to modify the code to use
c  other data intervals.
c 
c 
#endif
#include "def_master.h"
#include "param.h"
#include "cvbc.h"
c
      do 10 j=1,j_w
      do 10 i=1,i_w
        wstx(i,j) = wsxobs(j)
        wsty(i,j) = wsyobs(j)
  10  continue
c
      do 20 j=1,j_s
      do 20 i=1,i_s
        ssta(i,j) = sstobs(j)
        sssa(i,j) = salobs(j)
  20  continue
c
      return
      end
