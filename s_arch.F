      subroutine s_arch
c
c @(#) SCCS module: s_arch.F  version: 1.10
c     Creation date: 08/09/96
c
c==================================================================
c     Entry wrrest copies current model state into backup 
c     storage for asychronoustransfer to the master process.  
c     This transfer is done piecemeal in response to requests 
c     from the master.
c
c     Entry rdrest copies from the backup storage.
c
c     NOTE - this version only saves and restores the current
c     timestep.  The model must therefore be restarted with
c     a forward or euler backward timestep.
c==================================================================
c
#include "def_slave.h"
#include "param.h"
#include "ctmngr.h"
#include "frees.h"
#include "iounit.h"
#include "levind.h"
#include "slabs.h"
#include "switch.h"
#include "timelv.h"
#include "archive.h"
c
c------------------------------------------------------------------
c     read restart data set from archive for timesteps itt 
c------------------------------------------------------------------
c
      entry rdrest
c  check np equals 3
      if(np.ne.3)then
        write(ioslave,21) np, nc, nm
        stop
      endif
c
c------------------------------------------------------------------
c     initialize t,s,u,v.
c          velocities set to zero
c          temperature and salinity set to land values
c          -2.0 and 45.0 ppt
c------------------------------------------------------------------
c
c
      do 200 n=1,3
      do 200 j=1,JMT_S
      do 200 i=1,IMT_S
      do 200 k=1,KM
        u(k,i,j,n)=c0
        v(k,i,j,n)=c0
200   continue
      do 201 n=1,3
      do 201 j=1,JMT_S
      do 201 i=1,IMT_S
      do 201 k=1,KM
        t(k,i,j,1,n)= -2.0
        t(k,i,j,2,n)= 0.01
201   continue
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
c  transfer archived data from buffer storage.
c------------------------------------------------------------------
c
      if(mod(itt,nmix).ne.0)then
        write(ioslave,23) itt, nmix
        stop
      endif
c
      do 320 j=1,JMT_S
      do 320 i=1,IMT_S
       if(kmt(i,j).ne.0) then
         h0(i,j,3)= hsave(i,j,1)
       endif
       if(kmu(i,j).ne.0) then
         u0(i,j,3)= hsave(i,j,2)
         v0(i,j,3)= hsave(i,j,3)
       endif
320   continue
c
      do 350 n=1,NT
      do 350 j=1,JMT_S
      do 350 i=1,IMT_S
      do 350 k=1,KM
       if(k.le.kmt(i,j)) t(k,i,j,n,3)=tsave(k,i,j,n)
 350  continue
c
      do 380 j=1,JMT_S
      do 380 i=1,IMT_S
      do 380 k=1,KM
       if(k.le.kmu(i,j)) then
        u(k,i,j,3)=tsave(k,i,j,NT+1)
        v(k,i,j,3)=tsave(k,i,j,NT+2)
       endif
 380  continue
      return
c
  21  format(' Routine restio.',/,
     &       ' Trying to restart with np not equal to 3.',/,
     &       ' np, nc, nm  =',3i5,/,' Programme stopped.')
  23  format(' Routine restio - RESTART ERROR!',/,
     &       ' Trying to restart without forward timestep,',/,
     &       ' i.e. mod(itt,nmix).ne.0.    itt, nmix =',2i5,/,
     &       ' Programme stopped.')
c
c------------------------------------------------------------------
c     write restart data set to archive.
c------------------------------------------------------------------
c
      entry wrrest
c
c Copy model state information into assigned storage
c
      ittsav  = itt
      totssav = totsec
      totdsav = totday
      yrssav  = years
c
      do 400 j=1,JMT_S
      do 400 i=1,IMT_S
       hsave(i,j,1) = h0(i,j,np0)
       hsave(i,j,2) = u0(i,j,np0)
       hsave(i,j,3) = v0(i,j,np0)
400   continue
c
      do 450 n=1,NT
      do 450 j=1,JMT_S
      do 450 i=1,IMT_S
      do 450 k=1,KM
       tsave(k,i,j,n) = t(k,i,j,n,np)
 450  continue
c
      do 500 j=1,JMT_S
      do 500 i=1,IMT_S  
      do 500 k=1,KM 
       tsave(k,i,j,NT+1) = u(k,i,j,np)
       tsave(k,i,j,NT+2) = v(k,i,j,np)
 500  continue
      return
c
      end
