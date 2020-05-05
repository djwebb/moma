      subroutine ocn1st
c
c @(#) SCCS module: ocn1st.F  version: 1.1
c     Creation date: 03/16/95
c
c==================================================================
c
c  initialize quantites when model first started
c
c==================================================================
c
#include "def_slave.h"
#include "param.h"
#include "coord.h"
#include "ctmngr.h"
c#include "cvbc.h"
#include "frees.h"
#include "grdvar.h"
#include "iounit.h"
#include "levind.h"
#include "scalar.h"
#include "slabs.h"
#include "switch.h"
c
c------------------------------------------------------------------
c     generate topography over "t" points
c------------------------------------------------------------------
c
c     call topog !Topography now initialised in master program!
c
c------------------------------------------------------------------
c     initialize t,s,u,v.
c          velocities set to zero
c          temperature an salinity set to 0.0 and 35.0 at sea
c                                      to -2.0 and 45.0 over land.
c------------------------------------------------------------------
c
      sconst=0.0349
      do 200 n=1,3
      do 200 j=1,JMT_S
      do 200 i=1,IMT_S
      do 200 k=1,KM
        u(k,i,j,n)=c0
        v(k,i,j,n)=c0
200   continue
      do 201 k=1,KM
      do 201 j=1,JMT_S
      do 201 i=1,IMT_S
      do 201 n=1,3
       if(k.gt.kmt(i,j))then
           t(k,i,j,1,n)=-2.0
           t(k,i,j,2,n)=0.01
         else
           t(k,i,j,1,n)= theta0 (yt(j), zt(k))
           t(k,i,j,2,n)= sconst-0.035
       endif
201   continue
c
c------------------------------------------------------------------
c     initialize free surface model
c------------------------------------------------------------------
c
      do 300 n=1,3
      do 300 j=1,JMT_S
      do 300 i=1,IMT_S
        h0(i,j,n) = c0
        u0(i,j,n) = c0
        v0(i,j,n) = c0
300   continue
c
c------------------------------------------------------------------
c     initialize controls
c------------------------------------------------------------------
c
      itt   = 0
      totsec = c0
      totday = c0
      years  = c0
c
      return
      end
