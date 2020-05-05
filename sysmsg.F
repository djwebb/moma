      subroutine sysmsg
c
c @(#) SCCS module: sysmsg.F  version: 1.1
c     Creation date: 10/13/97
c
c==================================================================
c  standard pvm system routine calls for moma.pvm
c  subroutine is not used by the MPI system
c==================================================================
c
#include "pconst.h"
c     integer nids(nprc)
c     character slave*(*), arch*(*)
c
c------------------------------------------------------------------
c Spawn "nprc" copies of slave program
c------------------------------------------------------------------
c
c     entry spawn( slave, iflag, arch, nprc, nids, info )
c     call pvmfspawn(  slave, iflag, arch, nprc, nids, info )
c     return
c
c------------------------------------------------------------------
c Close down process and leave the system
c------------------------------------------------------------------
c
      entry closedown(info)
      if(me.eq.0)then
        write (stdout,14)
      endif
      call MPI_FINALIZE(ierr)
      stop
  14  format(1x,' ==> END of model run.',//)
c
c------------------------------------------------------------------
c Find the process id of this process
c------------------------------------------------------------------
c
c     entry getmyid(mytid)
c      call pvmfmytid(mytid)
c      if(mytid.lt.0)then
c        print *, ' pvmfmytid responds with negative value.',
c    &            ' program stopping'
c        stop
c      endif
c      return
c
c------------------------------------------------------------------
c Find out who spawned this process
c------------------------------------------------------------------
c
c      entry getmasterid(mtid)
c      call pvmfparent(mtid)
c      return
       end

