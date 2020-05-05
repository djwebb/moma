      subroutine sysmsg
c
c @(#) SCCS module: sysmsg.F  version: 1.1
c     Creation date: 03/16/95
c
c==================================================================
c  standard pvm system routine calls for moma.pvm
c==================================================================
c
#include "fpvm3.h"
      integer nids(nprc)
      character slave*(*), arch*(*)
c
c------------------------------------------------------------------
c Spawn "nprc" copies of slave program
c------------------------------------------------------------------
c
      entry spawn( slave, iflag, arch, nprc, nids, info )
      call pvmfspawn(  slave, iflag, arch, nprc, nids, info )
      return
c
c------------------------------------------------------------------
c Close down process and leave the system
c------------------------------------------------------------------
c
      entry closedown(info)
      call pvmfexit(info)
      return
c
c------------------------------------------------------------------
c Find the process id of this process
c------------------------------------------------------------------
c
      entry getmyid(mytid)
       call pvmfmytid(mytid)
       if(mytid.lt.0)then
         print *, ' pvmfmytid responds with negative value.',
     &            ' program stopping'
         stop
       endif
       return
c
c------------------------------------------------------------------
c Find out who spawned this process
c------------------------------------------------------------------
c
       entry getmasterid(mtid)
       call pvmfparent(mtid)
       return
       end

