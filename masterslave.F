     	program masterslave
c
c @(#) SCCS module: masterslave.F  version: 1.3
c     Creation date: 03/17/95
c
c==================================================================
c     top level routine for running moma.pvm on a multiprocessor
c     where the master and slave components have to be combined
c     into a single program.
c==================================================================
c
#include "fpvm3.h"
#include "def_master.h"
#include "param.h"
#include "mesdta.h"
#include "varinfo.h"
#include "ctmngr.h"
#include "snaps.h"
#include "levind.h"
#ifdef cray-t3d
      	intrinsic my_pe
        me=my_pe()
	time_start=rtc()
	utime(3)=time_start
#else
	integer time
	itime_start=time()
	utime(3)=dtime(utime)
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
		write(stdout,*)'calling master ', me
                write(stdout,*)'PVM_MAX_PACK = ', imaxpack
                call flush(stdout)
		call master
                call flush(stdout)
        else
		mtid=0
		write(stdout,*)'calling slave ', me
                call flush(stdout)
		call slave
                call flush(stdout)
        endif
	time_end=rtc()
	total_time=(time_end-time_start)*6.6667e-09
	print*,'total elapsed time for ',me,' is ',total_time
#else
	call getmasterid(mtid)
	if (mtid .lt. 0) then
		print*,'calling master'
		call master
	else
		print*,'calling slave'
		call slave
	endif
	itime_end=time()
	utime(3)=dtime(utime)+utime(3)
	itotal_time=itime_end-itime_start
	print*,'total elapsed time for ',me,' is ',itotal_time,utime(3)
#endif
        end

