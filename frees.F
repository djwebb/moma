#ifndef inline1
      subroutine frees(ic,jc)
#endif
#if defined inline1 || !defined inline
#ifndef inline1
c
c @(#) SCCS module: frees.F  version: 1.8
c      Creation date: 08/23/95
c
c====================================================================
c
c  Based on the routine in "Killworth, P.D., Stainforth, D.,
c  Webb, D.J. and Paterson, S.M. (1989).  A free surface
c  Bryan-Cox-Semtner model.  Institute of Oceanographic
c  Sciences, Report No. 270.  Wormley, Godlaming, U.K.. 184pp."
c
c  Later modified to use a leapfrog and time smoothing scheme,
c  suggested by Semtner, to prevent aliasing instabilities
c  near the equator
c
c=====================================================================
c
c---------------------------------------------------------------------
c  define global data
c---------------------------------------------------------------------
c
#include "def_slave.h"
#include "param.h"
c
#include "chmix.h"
#include "ctmngr.h"
#include "frees.h"
#include "grdvar.h"
#include "levind.h"
#include "scalar.h"
#include "switch.h"
#include "timelv.h"
#include "mesdta.h"
#else
      ic=ia(l)
      jc=ja(l)
#endif
c
c---------------------------------------------------------------------
c  set indices
c---------------------------------------------------------------------
c
      ip=ic+1
      im=ic-1
      jp=jc+1
      jm=jc-1
c
c---------------------------------------------------------------------
c  timestep height field.
c---------------------------------------------------------------------
c
      if(kmt(ic,jc).ne.0)then
        boxar = dxr*dyr*cstr(jc)
#ifdef de_checkbd
c
c Calculate masks for the del-cross and del-plus operators
c on selected time steps for free-surface calculation
c
       lchkbd = itbt.eq.2
       if (lchkbd) then
         m_N  = min(1,kmt(ic,jp))
         m_S  = min(1,kmt(ic,jm))
         m_E  = min(1,kmt(ip,jc))
         m_W  = min(1,kmt(im,jc))
         m_NE = min(1,kmu(ic,jc))
         m_SE = min(1,kmu(ic,jm))
         m_SW = min(1,kmu(im,jm))
         m_NW = min(1,kmu(im,jc))
c
c Calculate the del-plus operator (x2 for inclusion within existing terms
c of the dhdt expression)
c
         delplus = c2*(m_E*(h0(ip,jc,nm0) - h0(ic,jc,nm0))
     &                -m_W*(h0(ic,jc,nm0) - h0(im,jc,nm0))
     &                +m_N*(h0(ic,jp,nm0) - h0(ic,jc,nm0))
     &                -m_S*(h0(ic,jc,nm0) - h0(ic,jm,nm0)))
c
c Calculate the del-cross operator (x2 for inclusion within existing terms
c of the dhdt expression)
c
         delcross = ( m_NE*(h0(ip,jp,nm0) - h0(ic,jc,nm0))
     &               -m_SW*(h0(ic,jc,nm0) - h0(im,jm,nm0))
     &               +m_NW*(h0(im,jp,nm0) - h0(ic,jc,nm0))
     &               -m_SE*(h0(ic,jc,nm0) - h0(ip,jm,nm0)))
       endif
c
c Need a constant of the form: 
c        WGHT = alpha*grav*Hmax*dtbt/(dphi*dlambda*radius*radius*cos(phi)) 
c    i.e.WGHT = alpha*grav*Hmax*dtbt*boxar = dchkbd*boxar
c
# endif
        dhdt  = (dy*(u0(im,jm,nc0)*h(im,jm)+u0(im,jc,nc0)*h(im,jc)
     &              -u0(ic,jm,nc0)*h(ic,jm)-u0(ic,jc,nc0)*h(ic,jc))
     & +dx*(csu(jm)*(v0(im,jm,nc0)*h(im,jm)+v0(ic,jm,nc0)*h(ic,jm))
     &     -csu(jc)*(v0(im,jc,nc0)*h(im,jc)+v0(ic,jc,nc0)*h(ic,jc))
     &     ))*0.5*boxar
# ifdef de_checkbd
       if(lchkbd) then
        dhdt=dhdt+dchkbd*(delplus-delcross)*0.5*boxar
       endif
# endif
c
        h0(ic,jc,np0) = h0(ic,jc,nm0) + c2dtbt*dhdt
      endif
c
c---------------------------------------------------------------------
c  timestep velocity fields using semi-implicit scheme
c  for the coriolis term.
c---------------------------------------------------------------------
c
      if(kmu(ic,jc).ne.0)then
c
c---------------------------------------------------------------------
c    compute pressure gradients
c---------------------------------------------------------------------
c
        fxa = grav*dx2r*csur(jc)
        fxb = grav*dy2r
        temp1 = h0(ip,jp,nc0) - h0(ic,jc,nc0)
        temp2 = h0(ic,jp,nc0) - h0(ip,jc,nc0)
        dpdx0 = (temp1-temp2)*fxa
        dpdy0 = (temp1+temp2)*fxb
#ifdef free_eb
c
c---------------------------------------------------------------------
c  initial euler backward scheme treats coriolis term implicitly
c---------------------------------------------------------------------
c
        fac1=fcor(jc)*c2dtbt*p5
        fac2=c1/(c1+fac1*fac1)
        ustar = u0(ic,jc,nm0) + c2dtbt*(zu(ic,jc) - dpdx0)
     &                        + fac1*v0(ic,jc,nm0)
        vstar = v0(ic,jc,nm0) + c2dtbt*(zv(ic,jc) - dpdy0)
     &                        - fac1*u0(ic,jc,nm0)
        u0(ic,jc,np0) = (ustar + fac1*vstar)*fac2
        v0(ic,jc,np0) = (vstar - fac1*ustar)*fac2
#else
c
c---------------------------------------------------------------------
c  leapfrog scheme (explicit euler-backward on mixing timestep)
c---------------------------------------------------------------------
c
        u0(ic,jc,np0) = u0(ic,jc,nm0) + c2dtbt*(zu(ic,jc) - dpdx0
     &                        + fcor(jc)*v0(ic,jc,nc0) )
        v0(ic,jc,np0) = v0(ic,jc,nm0) + c2dtbt*(zv(ic,jc) - dpdy0
     &                        - fcor(jc)*u0(ic,jc,nc0) )
c
c---------------------------------------------------------------------
c  Calculate time averages (except for first pass of euler-backward)
c---------------------------------------------------------------------
c
        if(itbtp.ne.0)then
          freeav(1,ic,jc) = freeav(1,ic,jc) + h0(ic,jc,np0)
          freeav(2,ic,jc) = freeav(2,ic,jc) + u0(ic,jc,np0)
          freeav(3,ic,jc) = freeav(3,ic,jc) + v0(ic,jc,np0)
          if(itbt.eq.ntbt2)then
            fx = 0.5/ntbt
            h0(ic,jc,np0) = freeav(1,ic,jc)*fx
            u0(ic,jc,np0) = freeav(2,ic,jc)*fx
            v0(ic,jc,np0) = freeav(3,ic,jc)*fx
          endif
        endif
#endif
      endif
c
#endif
#ifndef inline1
      return
      end
#endif
