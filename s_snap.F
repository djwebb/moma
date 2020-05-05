      subroutine snap
c
c @(#) SCCS module: s_snap.F  version: 1.10
c     Creation date: 08/09/96
c
c==================================================================
c
c     Copy current model state into assigned storage for asychronous
c     transfer to the master process. This transfer is done piecemeal
c     in response to requests from the master.
c
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
#include "snaps.h"
#include "mesdta.h"
c
c Copy model state information into assigned storage
c
      ittsnap  = itt
c
      do 500 iisnp = 1, NSNAPS
      kvar=ispvar(iisnp)
      klev=isplev(iisnp)
c
      if(kvar.eq.1) then
c
        do 400 j=1,JMT_S
        do 400 i=1,IMT_S
         snap2d(i,j,iisnp) = h0(i,j,np0)
400     continue
c
      elseif(kvar.eq.2) then
c
        do 410 j=1,JMT_S
        do 410 i=1,IMT_S
         snap2d(i,j,iisnp) = u0(i,j,np0)
410     continue
c
      elseif(kvar.eq.3) then
c
        do 420 j=1,JMT_S
        do 420 i=1,IMT_S
         snap2d(i,j,iisnp) = v0(i,j,np0)
420     continue
c
      elseif(kvar.eq.4) then
c
       do 430 j=1,JMT_S
       do 430 i=1,IMT_S
        snap2d(i,j,iisnp) = t(klev,i,j,1,np)
 430   continue
c
      elseif(kvar.eq.5) then
c
       do 440 j=1,JMT_S
       do 440 i=1,IMT_S
        snap2d(i,j,iisnp) = t(klev,i,j,2,np)*1000. + 35.0
 440   continue
c
      elseif(kvar.eq.6) then
c
       do 450 j=1,JMT_S
       do 450 i=1,IMT_S
        snap2d(i,j,iisnp) = u(klev,i,j,np)
 450   continue
c
      elseif(kvar.eq.7) then
c
       do 460 j=1,JMT_S
       do 460 i=1,IMT_S
        snap2d(i,j,iisnp) = v(klev,i,j,np)
 460   continue
c
      endif
c
 500  continue
      return
      end
