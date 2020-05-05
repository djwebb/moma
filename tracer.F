      subroutine tracer(ic,jc)
c
c @(#) SCCS module: tracer.F  version: 1.2
c     Creation date: 03/17/95
c
c=============================================================
c
c     tracer computes tracer quantities for time level t+1
c
c=============================================================
c
#include "def_slave.h"
#include "param.h"
c
#include "scalar.h"
#include "grdvar.h"
#include "levind.h"
#include "timelv.h"
#include "slabs.h"
c
#include "chmix.h"
c#include "cvbc.h"
#include "cvmix.h"
c	TASK COMMON /cvbc/ smf(2),bmf(2),stf(nt),btf(nt)
	COMMON /cvbc/ smf(2),bmf(2),stf(nt),btf(nt)
c	dimension fue(km),fuw(km),fvn(km),fvs(km),
c    &		  fw(0:km),fwb1(0:km),temp(km),
c    &		  dpdx(km),dpdy(km),vtf(0:km),
c    & 		  maskpo(km),maskmo(km),maskop(km),maskom(km)
c
c=============================================================
c    outer loops over i and j indices
c=============================================================
c
      kmc= kmt(ic,jc)
      if(kmc.eq.0)return
      jp = jc+1
      jm = jc-1
      ip = ic+1
      im = ic-1
      boxa = dx*cst(jc)*dy
      boxar = c1/boxa
c
c=============================================================
c    set up masks, equal to 0 on sub-seabottom points.
c=============================================================
c
      do 110 k=1,kmc
      maskop(k)=0
      if(k.le.kmt(ic,jp))maskop(k)=1
      maskom(k)=0
      if(k.le.kmt(ic,jm))maskom(k)=1
      maskpo(k)=0
      if(k.le.kmt(ip,jc))maskpo(k)=1
      maskmo(k)=0
      if(k.le.kmt(im,jc))maskmo(k)=1
 110  continue 
c
c=============================================================
c  calculate horizontal and vertical advective fluxes
c=============================================================
c
      dxp5 = p5*dx
      dyp5 = p5*dy
      do 160 k=1,kmc
      fvn(k) = (v(k,ic,jc,nc)+v(k,im,jc,nc))*dxp5*csu(jc)*dz(k)
      fvs(k) = (v(k,ic,jm,nc)+v(k,im,jm,nc))*dxp5*csu(jm)*dz(k)
      fue(k) = (u(k,ic,jc,nc)+u(k,ic,jm,nc))*dyp5*dz(k)
      fuw(k) = (u(k,im,jc,nc)+u(k,im,jm,nc))*dyp5*dz(k)
      fw(k-1)  = (fuw(k)-fue(k)+fvs(k)-fvn(k))
 160  continue

c
c=============================================================
c    integrate up from the lowest layer 
c    fw is positive upwards
c=============================================================
c
      fw(kmc)=c0
      do 180 k=kmc-1,0,-1
      fw(k)=fw(k)+fw(k+1)
 180  continue
c
c=============================================================
c    main calculation loop to timestep tracers
c=============================================================
c
      do 900 m=1,nt
c
c=============================================================
c  calculate vertical fluxes - positive upwards
c=============================================================
c
      do 930 k=1,kmc-1
        vtf(k) = fkph*(t(k+1,ic,jc,m,nm)-t(k  ,ic,jc,m,nm))
     &                 *dzwr(k)*boxa
        fwb1(k) = fw(k)*(t(k,ic,jc,m,nc)+t(k+1,ic,jc,m,nc))*p5
 930  continue
c
c=============================================================
c  surface and bottom advective and diffusive boundary conditions
c  stf is downwards flux into the ocean
c=============================================================
c
      vtf(0)   = -stf(m)*boxa
      vtf(kmc) = c0
      fwb1(0)  = fw(0)*t(1,ic,jc,m,nc)
      fwb1(kmc)= c0
c
c=============================================================
c  main timestepping calculation
c=============================================================
c
c
      do 960 k=1,kmc
      t(k,ic,jc,m,np) = t(k,ic,jc,m,nm)
     &                + c2dtts*boxar*dzr(k)*(
c  u-advection
     &   fuw(k)*(t(k,im,jc,m,nc) + t(k,ic,jc,m,nc))*p5
     & - fue(k)*(t(k,ip,jc,m,nc) + t(k,ic,jc,m,nc))*p5
c  v-advection
     & + fvs(k)*(t(k,ic,jm,m,nc) + t(k,ic,jc,m,nc))*p5
     & - fvn(k)*(t(k,ic,jp,m,nc) + t(k,ic,jc,m,nc))*p5
c  w-advection
     & + fwb1(k  ) - fwb1(k-1)
c  x diffusion
     & +(bbt(jc)*(t(k,ip,jc,m,nm)-t(k,ic,jc,m,nm))*maskpo(k)
     & + bbt(jc)*(t(k,im,jc,m,nm)-t(k,ic,jc,m,nm))*maskmo(k)
c  y diffusion
     & + cct(jc)*(t(k,ic,jp,m,nm)-t(k,ic,jc,m,nm))*maskop(k)
     & + ddt(jc)*(t(k,ic,jm,m,nm)-t(k,ic,jc,m,nm))*maskom(k))
     & * dz(k)
c  z diffusion  
     & + vtf(k  ) - vtf(k-1)
     &   )
 960  continue
 900  continue
c
c=============================================================
c   convectively adjust water column if gravitaionally unstable.
c=============================================================
c
      if(kmc.gt.1.and.ncon.gt.0)then
        do 1070 n=1,ncon
        do 1070 ks=1,2
        call state(t(1,ic,jc,1,np),t(1,ic,jc,2,np),temp(1),
     &              kmc, ks)
        do 1070 m=1,nt
        do 1070 k=ks,kmc-1,2
        if(temp(k).gt.temp(k+1))then
          t(k,ic,jc,m,np) = (dz(k)*t(k,ic,jc,m,np)
     &                    +  dz(k+1)*t(k+1,ic,jc,m,np))*dzw2r(k)
          t(k+1,ic,jc,m,np) = t(k,ic,jc,m,np)
        endif
 1070   continue
      endif
c
      return
      end
