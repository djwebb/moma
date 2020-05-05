     	
      program masterslave
c
c @(#) SCCS module: masterslave.F  version: 1.1
c     Creation date: 10/13/97
c
c==================================================================
c     top level routine for running moma.mpi on a multiprocessor
c     where the master and slave components have to be combined
c     into a single program. 
c==================================================================
c
#include "def_master.h"
c#include "mpif.h"
#include "param.h"
#include "mesdta.h"
#include "varinfo.h"
#include "ctmngr.h"
#include "snaps.h"
#include "levind.h"
c
c  use MPI timing routine for times
c
      double precision total_time
c
c  initialise MPI, initialise timing and request process number
c
      call MPI_INIT(ierr)
      start_time = MPI_WTIME()
      call MPI_BUFFER_ATTACH(bufmpi,LENBUF*NUMBUF,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,me,ierr)
c
c  call master or slave routine
c
      if (me .eq. 0) then
        write(stdout,10)me
        call master
      else
        write(stdout,11)me
        call slave
      endif
c
c  write out total time and stop
c
      total_time = MPI_WTIME() - start_time
      write(stdout,12)me,total_time
      if(me.eq.0)then
        write (stdout,14)me
      endif
c
      call MPI_FINALIZE(ierr)
      stop
c
  10  format('processor ',i4,' calling subroutine master')
  11  format('processor ',i4,' calling subroutine slave')
  12  format('total elapsed time for ',i4,' is ',f10.4)
  14  format(1x,' ==> END of model run. [masterslave exit]',i5,//)
c
      end


