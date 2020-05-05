      subroutine clinic(ic,jc)
c
c @(#) SCCS module: clinic.F  version: 1.2
c     Creation date: 03/17/95
c
c==================================================================
c
c     clinic computes internal mode velocities at time level t+1
c     for a vertical column of k cells with i and j values equal to
c     ic and jc.
c
c==================================================================
c
c
#include "def_slave.h"
#include "param.h"
c
#include "scalar.h"
#include "switch.h"
#include "coord.h"
#include "grdvar.h"
#include "levind.h"
#include "timelv.h"
#include "slabs.h"
#include "frees.h"
#include "mesdta.h"
c
#include "chmix.h"
c#include "cvbc.h"
#include "cvmix.h"
c
c	TASK COMMON /cvbc/ smf(2),bmf(2),stf(nt),btf(nt)
	COMMON /cvbc/ smf(2),bmf(2),stf(nt),btf(nt)
c	dimension fue(km),fuw(km),fvn(km),fvs(km),
c    &		  fw(0:km),fwb1(0:km),fwb2(0:km),
c    &		  dpdx(km),dpdy(km),vmf(0:km,2),
c    &		  rhoo(km),rhpo(km),rhpp(km),rhop(km)

      kmc=kmu(ic,jc)
      if(kmc.eq.0)return
      jp = jc+1
      jm = jc-1
      ip = ic+1
      im = ic-1
c
c  find max of all the neighbouring points.
c
c#ifdef oldadv
c      kmd  = max(kmu(ip,jc),kmu(im,jc),kmu(ic,jp),kmu(ic,jm),kmc)
c#else
c      kmd  = max(kmu(im,jm),kmu(ic,jm),kmu(ip,jm),
c     &           kmu(im,jc),kmc       ,kmu(ip,jc),
c     &           kmu(im,jp),kmu(ic,jp),kmu(ip,jp))
c#endif
      boxa = dx*csu(jc)*dy
      boxar = c1/boxa
c
c=============================================================
c  calculate horizontal and vertical advective fluxes
c=============================================================
c
      dxp5 = p5*dx
      dyp5 = p5*dy
      do 100 k=1,kmd(ic,jc)
#ifdef oldadv
c
c  original horizontall advection scheme
c
      fvn(k) = (v(k,ic,jc,nc)+v(k,ic,jp,nc))*p5*dx*cst(jp)*dz(k)
      fvs(k) = (v(k,ic,jc,nc)+v(k,ic,jm,nc))*p5*dx*cst(jc)*dz(k)
      fue(k) = (u(k,ic,jc,nc)+u(k,ip,jc,nc))*p5*dy*dz(k)
      fuw(k) = (u(k,ic,jc,nc)+u(k,im,jc,nc))*p5*dy*dz(k)
#else
c
c  new horizontal advection scheme
c
      fvjc   = csu(jc)*
     &         (v(k,im,jc,nc)+2.0*v(k,ic,jc,nc)+v(k,ip,jc,nc))
      fuic   = (u(k,ic,jp,nc)+2.0*u(k,ic,jc,nc)+u(k,ic,jm,nc))
c
      fvn(k) = dz(k)*dx*p125*( fvjc + csu(jp)*
     &         (v(k,im,jp,nc)+2.0*v(k,ic,jp,nc)+v(k,ip,jp,nc)))
      fvs(k) = dz(k)*dx*p125*( fvjc + csu(jm)*
     &         (v(k,im,jm,nc)+2.0*v(k,ic,jm,nc)+v(k,ip,jm,nc)))
c
      fue(k) = dz(k)*dy*p125*( fuic +
     &         (u(k,ip,jp,nc)+2.0*u(k,ip,jc,nc)+u(k,ip,jm,nc)))
      fuw(k) = dz(k)*dy*p125*( fuic +
     &         (u(k,im,jp,nc)+2.0*u(k,im,jc,nc)+u(k,im,jm,nc)))
#endif
c
c  new vertical boundary conditions on w at u points
c
      fw(k-1)  = (fuw(k)-fue(k)+fvs(k)-fvn(k))
 100  continue
c
c=============================================================
c    integrate up from the bottom level
c=============================================================
c
      fw(kmd(ic,jc))=c0
      do 210 k=kmd(ic,jc),1,-1
      fw(k-1)=fw(k-1)+fw(k)
 210  continue
#ifdef presetp
c
c-----------------------------------------------------------------
c         compute hydrostatic pressure gradient (baroclinic part)
c-----------------------------------------------------------------
c
      fxa = dx2r*csur(jc)
      fxb = dy2r
      do 300 k=1,kmc
        temp1  = p(k,ip,jp) - p(k,ic,jc)
        temp2  = p(k,ic,jp) - p(k,ip,jc)
        dpdx(k) = (temp1-temp2)*fxa
        dpdy(k) = (temp1+temp2)*fxb
 300  continue
#else
c
c-----------------------------------------------------------------
c     compute density on surrounding 't' points
c-----------------------------------------------------------------
c
      call state (t(1,ic,jc,1,nc), t(1,ic,jc,2,nc), rhoo, 
     &            kmt(ic,jc) ,0)
      call state (t(1,ip,jc,1,nc), t(1,ip,jc,2,nc), rhpo, 
     &            kmt(ip,jc) ,0)
      call state (t(1,ip,jp,1,nc), t(1,ip,jp,2,nc), rhpp, 
     &            kmt(ip,jp) ,0)
      call state (t(1,ic,jp,1,nc), t(1,ic,jp,2,nc), rhop, 
     &            kmt(ic,jp) ,0)
c
c-----------------------------------------------------------------
c         compute hydrostatic pressure gradient (baroclinic part)
c-----------------------------------------------------------------
c
c    1. compute it at the first level
c
      fxa = grav*dzw(0)*dx2r*csur(jc)
      fxb = grav*dzw(0)*dy2r
      temp1  = rhpp(1) - rhoo(1)
      temp2  = rhop(1) - rhpo(1)
      dpdx(1) = (temp1-temp2)*fxa
      dpdy(1) = (temp1+temp2)*fxb
c
c    2. compute the change in pressure gradient between levels
c
      fxa = grav*dx4r*csur(jc)
      fxb = grav*dy4r
      do 300 k=2,kmc
        temp1 = (rhpp(k)+rhpp(k-1))- (rhoo(k)+rhoo(k-1))
        temp2 = (rhop(k)+rhop(k-1))- (rhpo(k)+rhpo(k-1))
        dpdx(k) = (temp1-temp2)*fxa*dzw(k-1)
        dpdy(k) = (temp1+temp2)*fxb*dzw(k-1)
 300  continue
c
c=============================================================
c         integrate downward from the first level
c=============================================================
c
      do 400 k=2,kmc
        dpdx(k) = dpdx(k-1) + dpdx(k)
        dpdy(k) = dpdy(k-1) + dpdy(k)
 400  continue
#endif
c
c-----------------------------------------------------------------
c         calculate quantities for the computation of
c         vertical diffusion and advection of momentum
c         - positive for upwards flux of momentum
c-----------------------------------------------------------------
c
      do 500 k=1,kmc-1
        vmf(k,1) = fkpm*(u(k+1,ic,jc,nm) - u(k  ,ic,jc,nm))
     &                     *dzwr(k)*dx*csu(jc)*dy
        vmf(k,2) = fkpm*(v(k+1,ic,jc,nm) - v(k  ,ic,jc,nm))
     &                     *dzwr(k)*dx*csu(jc)*dy
        fwb1(k) = fw(k)*(u(k,ic,jc,nc) + u(k+1,ic,jc,nc))*p5
        fwb2(k) = fw(k)*(v(k,ic,jc,nc) + v(k+1,ic,jc,nc))*p5
 500  continue
c
c-----------------------------------------------------------------
c     set k=0 vertical diffusion to reflect surface wind stress
c     set level "kmc" for bottom drag condition.
c     set k=kmc elements of vertical advection array to zero
c     vertical momentum transfer.  Set k=0 to represent advection
c     loss to region above the free surface.
c     smf is the surface stress acting on the ocean
c     bmf is the bottom stress acting on the ocean floor
c-----------------------------------------------------------------
c
      vmf(0,1)  = -smf(1)*boxa
      vmf(0,2)  = -smf(2)*boxa
      vmf(kmc,1)= -bmf(1)*boxa
      vmf(kmc,2)= -bmf(2)*boxa
      fwb1(0)   = fw(0)*u(1,ic,jc,nc)
      fwb2(0)   = fw(0)*v(1,ic,jc,nc)
      fwb1(kmc) = c0
      fwb2(kmc) = c0
c
c------------------------------------------------------------------
c       set up index lcor for the coriolis term:
c       use tau timestep for explicit trtmnt ... or
c       use tau-1 timestep for implicit treatment with remainder of
c       term to be added later
c------------------------------------------------------------------
c
          if (acor .eq. c0) then
            lcor = nc
          else
            lcor = nm
          endif
c
c------------------------------------------------------------------
c       initialise vertical mean forcing array elements:
c------------------------------------------------------------------
c
      zu(ic,jc) = c0
      zv(ic,jc) = c0
c
c=============================================================
c    main calculation to timestep baroclinic velocity
c=============================================================
c
      do 600 k=1,kmc
c
c u-momentum equation.  separate coriolis term and the rest.
c
      ucor  = fcor(jc)*v(k,ic,jc,lcor)
      urest = -dpdx(k)
c  u-advection
     & +(fuw(k)*(u(k,im,jc,nc) + u(k,ic,jc,nc))*p5
     & - fue(k)*(u(k,ip,jc,nc) + u(k,ic,jc,nc))*p5
c  v-advection
     & + fvs(k)*(u(k,ic,jm,nc) + u(k,ic,jc,nc))*p5
     & - fvn(k)*(u(k,ic,jp,nc) + u(k,ic,jc,nc))*p5
c  w-advection
     & + fwb1(k  ) - fwb1(k-1)
c  x diffusion  
     & +(bbu(jc)*(u(k,ip,jc,nm)-u(k,ic,jc,nm))
     & + bbu(jc)*(u(k,im,jc,nm)-u(k,ic,jc,nm))
c  y diffusion
     & + ccu(jc)*(u(k,ic,jp,nm)-u(k,ic,jc,nm))
     & + ddu(jc)*(u(k,ic,jm,nm)-u(k,ic,jc,nm)) )
     & *dz(k)
c  z diffusion  
     &  + vmf(k  ,1) - vmf(k-1,1) )*boxar*dzr(k)
c
c v-momentum equation.  separate coriolis term and the rest
c
      vcor  = -fcor(jc)*u(k,ic,jc,lcor)
      vrest = -dpdy(k)
c  u-advection
     & +(fuw(k)*(v(k,im,jc,nc) + v(k,ic,jc,nc))*p5
     & - fue(k)*(v(k,ip,jc,nc) + v(k,ic,jc,nc))*p5
c  v-advection
     & + fvs(k)*(v(k,ic,jm,nc) + v(k,ic,jc,nc))*p5
     & - fvn(k)*(v(k,ic,jp,nc) + v(k,ic,jc,nc))*p5
c  w-advection
     & + fwb2(k  ) - fwb2(k-1)
c  x diffusion  
     & +(bbu(jc)*(v(k,ip,jc,nm)-v(k,ic,jc,nm))
     & + bbu(jc)*(v(k,im,jc,nm)-v(k,ic,jc,nm))
c  y diffusion
     & + ccu(jc)*(v(k,ic,jp,nm)-v(k,ic,jc,nm))
     & + ddu(jc)*(v(k,ic,jm,nm)-v(k,ic,jc,nm)) )
     & *dz(k)
c  z diffusion
     & + vmf(k  ,2) - vmf(k-1,2) )*boxar*dzr(k)
c
c***************************************************************
c  the corrections to the advection and diffusion terms that
c  allow for the curvature of the grid near the poles have been
c  left out of this version of the model
c***************************************************************
c
c-----------------------------------------------------------------
c     calculate vertically averaged forcing (without coriolis)
c     and total force at each level (with coriolis)
c-----------------------------------------------------------------
c
      fx = dz(k)
      zu(ic,jc) = zu(ic,jc) + urest*fx
      zv(ic,jc) = zv(ic,jc) + vrest*fx
      u(k,ic,jc,np) = ucor + urest
      v(k,ic,jc,np) = vcor + vrest
 600  continue
c
c-----------------------------------------------------------------
c     calculate vertical mean force 
c-----------------------------------------------------------------
c
      zu(ic,jc) = zu(ic,jc)*hr(ic,jc)
      zv(ic,jc) = zv(ic,jc)*hr(ic,jc)
c
c-----------------------------------------------------------------
c     timestep baroclinic terms
c-----------------------------------------------------------------
c
      if(acor.eq.c0)then
        do 700 k=1,kmc
          u(k,ic,jc,np) = u(k,ic,jc,nm) + c2dtuv*u(k,ic,jc,np) 
          v(k,ic,jc,np) = v(k,ic,jc,nm) + c2dtuv*v(k,ic,jc,np) 
 700    continue
      else
        fx    = c2dtuv*acor*c2*omega*sine(jc)
        detmr = c1/(c1 + fx*fx)
        do 750 k=1,kmc
          t1 = (u(k,ic,jc,np)+fx*v(k,ic,jc,np))*detmr
          t2 = (v(k,ic,jc,np)-fx*u(k,ic,jc,np))*detmr
          u(k,ic,jc,np) = u(k,ic,jc,nm) + c2dtuv*t1 
          v(k,ic,jc,np) = v(k,ic,jc,nm) + c2dtuv*t2
 750    continue
      endif
c
c-----------------------------------------------------------------
c     remove vertical mean velocity
c-----------------------------------------------------------------
c
      sumu = c0
      sumv = c0
      zwr = c1/zw(kmc)
      do 800 k=1,kmc
      sumu = sumu + u(k,ic,jc,np)*dz(k)
      sumv = sumv + v(k,ic,jc,np)*dz(k)
 800  continue
      sumu = sumu*zwr
      sumv = sumv*zwr
      do 900 k=1,kmc
      u(k,ic,jc,np) = u(k,ic,jc,np) - sumu
      v(k,ic,jc,np) = v(k,ic,jc,np) - sumv
 900  continue
      return
      end
