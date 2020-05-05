     	program masterslave
c
c @(#) SCCS module: masterslave.F  version: 1.9
c     Creation date: 12/18/95
c
c==================================================================
c     top level routine for running moma.pvm on a multiprocessor
c     where the master and slave components have to be combined
c     into a single program.
c==================================================================
c
#include "def_master.h"
#include "param.h"
#include "mesdta.h"
#include "varinfo.h"
#include "ctmngr.h"
#include "snaps.h"
#include "levind.h"
#ifdef cray-t3d
      	intrinsic my_pe
        me = my_pe()
	time_start=rtc()*6.6e-9
        intime = time_start
#else
	integer time
	intime = time()
#endif	
	call getmyid(mytid)
#ifdef cray-t3d
c
c get the PVM_MAX_PACK environment setting and convert to words
c
        call pvmfgetopt(PVMMAXPACK,imaxpack)
        imaxpack=imaxpack/NBYTER
c
	do 10 i=1,n$pes
		itida(i)=i
 10	continue
	mytid=me

        if (me .eq. 0) then
	  write(stdout,*)'processor ', me,' calling subroutine master'
          write(stdout,*)'PVM_MAX_PACK = ', imaxpack
          call flush(stdout)
	  call master
          call flush(stdout)
        else
          mtid=0
	  write(stdout,*)'processor ', me,' calling subroutine slave'
          call flush(stdout)
	  call slave
          call flush(stdout)
        endif
        total_time=rtc()*6.6e-9-time_start
	print*,'total elapsed time for ',me,' is ',total_time
#else
	call getmasterid(mtid)
	if (mtid .lt. 0) then
	  write(stdout,*)'processor ', me,' calling subroutine master'
	  call master
	else
	  write(stdout,*)'processor ', me,' calling subroutine slave'
	  call slave
	endif
	utime(3) = etime(utime)
	itotal_time = time() - intime
	print * ,'total elapsed time for ',me,' is ',
     &            itotal_time,utime(3)
#endif
        end

